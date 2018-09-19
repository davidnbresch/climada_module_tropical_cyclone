function result = calibrate_TC_DF_emdat_region(x,fixed_parameters,delta_shape_parameter,entity,hazard,em_data,norm,bounds,type)
% NAME: calibrate_TC_DF_emdat_region
% MODULE: tropical_cyclone
%
%   This function (1) changes the TC damagefunction based on v_threshold =x(1), delta_shape and scale=x(2), 
%   (2) calculates the YDS based on entity, hazard and the new damage
%   function and (3) provides the R-squared difference between the YDS and
%   em_data (or other metric if requested)
% CALLING SEQUENCE:
%   Called as an objective function from function calibrate_TC_DF_emdat_framweork_region. See this
%   function for more info.
% EXAMPLE:
%
%   x0 = [30 .5]; % starting values of free parameters (v_threshold and scale)
%   delta_shape_parameter = 49; % difference between v_threshold and v_half
%   bounds.lb = [25 1e-9]; % never set lower bound of scale to 0! Otherwise calibration goes wrong.
%   bounds.ub = [35 1.];
%   norm.lb = ones(size(bounds.lb));
%   norm.ub = norm.lb+1;
%   norm.x0 = (x0-bounds.lb)./(bounds.ub-bounds.lb) .* (norm.ub-norm.lb) + norm.lb;
% 
%   fun = @(x)calibrate_TC_DF_emdat_region(x,delta_shape_parameter, entity, hazard, em_data_region, norm, bounds);
%   
%   options = optimoptions('patternsearch','UseParallel',false,...
%     'UseCompletePoll', true, 'UseVectorized', false,...
%     'MaxFunctionEvaluations',1200,'Display','iter',...
%     'Cache','on','InitialMeshSize',.25,...
%     'PollMethod','GPSPositiveBasis2N','StepTolerance',0.001);
%     
%   [x_result,fval] = patternsearch(fun,norm.x0,[],[],[],[],norm.lb,norm.ub,[],options);
%
% INPUTS:
%   x: double/vector, scale parameter (and optionally) shape_param to be
%   calibrated (normalized or not): [v_threshold scale]
%   !!! If norm and bounds are provided, x is assumed to be normalized.
%   
% OPTIONAL INPUTS:
%   delta_shape_parameter: difference im m/s between v_half and v_threshold
%   entity: struct or file name
%   hazard: struct
%   em_data: struct with fields year and damage
%   norm: struct normalized bounds and starting value x0 
%   bounds: struct with fields lb and ub, 2 values each
%   type (string):
%           'R2':   (DEFAULT) result is R^2 (= the sum of squared differences of
%                   year damages of emdat and climada for each specific historical year).
%           'AED':  "Annual Expected Damage": result is the squared
%                   difference of mean year damage 
%           'R':    result is R (= the sum of the absolute differences of
%                   year damages of emdat and climada for each specific historical year).
%           'RP':   "Return Period": as AED but for different return
%                   periods with weights (not implemented yet) - only makes
%                   sense for long time series
%
% OUTPUT:
%   result: metric requested by input vraiable type. (by DEFAULT, result is R^2 = the sum of squared differences of
%                   year damages of emdat and climada for each specific historical year)
%
% MODIFICATION HISTORY:
% Samuel Eberenz, eberenz@posteo.eu, 20180718, init
% Samuel Eberenz, eberenz@posteo.eu, 20180905, handle entity file name instead of struct
% Samuel Eberenz, eberenz@posteo.eu, 20180914, add option of 1 free parameter
% Samuel Eberenz, eberenz@posteo.eu, 20180919, add optimizerType 'logR2'
%-


global climada_global
damage_at_centroid_temp = climada_global.damage_at_centroid;
climada_global.damage_at_centroid = 1;

if ~exist('x','var'),error('Not enough input parameters');end %
if ~exist('fixed_parameters','var'),fixed_parameters=[];end %
if ~exist('delta_shape_parameter','var'),delta_shape_parameter=49;end % m/s
if xor(exist('bounds','var'), exist('norm','var')) % check if variables for normalizatiuon are given as input, if no, assume non-normalized input x
    error('Missing input: Function requires either both optional input variables norm and bounds or neither.')
elseif (exist('norm','var')&& exist('bounds','var')) % compute x from normalized x using norm and bound:
    norm.x = x; % if x was given normalized, do inverse normalization to get actual value of x
    x=(norm.x-norm.lb).*(bounds.ub-bounds.lb)./(norm.ub-norm.lb)+bounds.lb;
end
% disp(x)
if ~exist('type','var'),type='R2';end 
if ~exist('entity','var'),error('Not enough input parameters');end 
if isempty(entity),error('variable entity not specified. provide struct or path');end
if ~isstruct(entity)
    entity = climada_entity_load(entity);
end
    
result = 0; %init

switch length(x)
    case 1
        if isempty(fixed_parameters)
            error('x and fixed_parameters mismatch');
        end
        v_threshold = fixed_parameters;
        scale = x;
        v_half = v_threshold + delta_shape_parameter; % v_half
    case 2 % v_threshold fixed
        v_threshold = x(1); % v_threshold
        scale = x(end);
        v_half = v_threshold + delta_shape_parameter; % v_half
      
    otherwise
        error('length(x) has to be 2 currently');
end
disp([v_threshold v_half scale])
if scale<0 || scale>100 || v_half<=v_threshold || v_threshold<0
        error('Error: scale or shape parameter out of bounds');
end

%%%% Damage Function after Emanuel (2011) (c.f. also Elliott et al. 2015, similar also Sealy & Strobl et al. 2017):
%% change damage function of entity according to scale and shape parameters:
% so far damage function TC 001 is changed (to be refined):
choose_DF = strcmp(entity.damagefunctions.peril_ID,'TC') .* (entity.damagefunctions.DamageFunID==1);

% set PPA to 1 for all intensities
entity.damagefunctions.PAA(choose_DF==1) = 1;
% set MDD with x(1)= v_half [m/s] % 
v_temporary = max((entity.damagefunctions.Intensity(choose_DF==1) - v_threshold),0)/(v_half-v_threshold);
entity.damagefunctions.MDD(choose_DF==1)= v_temporary.^3 ./ (1+v_temporary.^3);
% linearly dampen MDD by multiplication with scale x(2), 0<x(2)<=1
entity.damagefunctions.MDD(choose_DF==1)= scale * entity.damagefunctions.MDD(choose_DF==1);
clear v_temporary
%figure;climada_damagefunctions_plot(entity,'TC');
%% (2) calculate YDS
EDS = climada_EDS_calc(entity,hazard,[],[],2);

clear entity choose_DF x bounds norm VV_temp 
%% (3) provide difference to em_data
switch type
    case 'AED'
        em_data_yyyy_allYears = em_data.year(1):em_data.year(end);
        em_data_damage_allYears = zeros(size(em_data_yyyy_allYears));
        for year_i = em_data_yyyy_allYears
            if max(ismember(em_data.year,year_i))
                em_data_damage_allYears(em_data_yyyy_allYears==year_i)= em_data.damage(em_data.year==year_i);
            end
        end
        result = (mean(em_data_damage_allYears)-EDS.ED)^2;
    case 'R2'
        [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
        [LIA,LOCB] = ismember(em_data.year,YDS.yyyy);
        % exclude years with 0-damage-entries:
        year_i_vector = find(LIA);
        year_i_vector(em_data.damage(year_i_vector)==0)=[]; 
        result = sum((em_data.damage(year_i_vector) - YDS.damage(LOCB(year_i_vector))).^2);

    case 'logR2'
        [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
        [LIA,LOCB] = ismember(em_data.year,YDS.yyyy);
        % exclude years with 0-damage-entries:
        year_i_vector = find(LIA);
        year_i_vector(em_data.damage(year_i_vector)==0)=[]; 
        result = sum(log((em_data.damage(year_i_vector) - YDS.damage(LOCB(year_i_vector))).^2));
        
    case 'R4'
        [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
        [LIA,LOCB] = ismember(em_data.year,YDS.yyyy);
        % exclude years with 0-damage-entries:
        year_i_vector = find(LIA);
        year_i_vector(em_data.damage(year_i_vector)==0)=[]; 
      %  year_i_vector(YDS.damage(LOCB(year_i_vector))==0)=[];

        result = sum((em_data.damage(year_i_vector) - YDS.damage(LOCB(year_i_vector))).^4);
    case 'R'
                [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
        [LIA,LOCB] = ismember(em_data.year,YDS.yyyy);
        % exclude years with 0-damage-entries:
        year_i_vector = find(LIA);
        year_i_vector(em_data.damage(year_i_vector)==0)=[]; 
    %    year_i_vector(YDS.damage(LOCB(year_i_vector))==0)=[];

        result = sum(abs(em_data.damage(year_i_vector) - YDS.damage(LOCB(year_i_vector))));
        
    case 'logR'
                [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
        [LIA,LOCB] = ismember(em_data.year,YDS.yyyy);
        % exclude years with 0-damage-entries:
        year_i_vector = find(LIA);
        year_i_vector(em_data.damage(year_i_vector)==0)=[]; 
     %   year_i_vector(YDS.damage(LOCB(year_i_vector))==0)=[];

        result = sum(log(abs(em_data.damage(year_i_vector) - YDS.damage(LOCB(year_i_vector)))));

    case 'RP'
        warning('Type RP not yet implemented, returning squared difference in AED instead');
        em_data_yyyy_allYears = em_data.year(1):em_data.year(end);
        em_data_damage_allYears = zeros(size(em_data_yyyy_allYears));
        for year_i = em_data_yyyy_allYears
            if max(ismember(em_data.year,year_i))
                em_data_damage_allYears(em_data_yyyy_allYears==year_i)= em_data.damage(em_data.year==year_i);
            end
        end
        result = (mean(em_data_damage_allYears)-EDS.ED)^2;
    otherwise
        error(['Unidentified optimizer type: ' typr]);
end

climada_global.damage_at_centroid = damage_at_centroid_temp;
clear EDS YDS em_data LIA LOCB year_i damage_at_centroid_temp

end

