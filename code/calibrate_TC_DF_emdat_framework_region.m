function result=calibrate_TC_DF_emdat_framework_region(TCBasinID,value_mode,cropped_assets,resolution,calibrate_countries,hazard_filename,number_free_parameters,years_considered,on_cluster,hand_over_entity_file)
% NAME: calibrate_TC_DF_emdat_framework_region
% MODULE: tropical_cyclone
%
%   This function calibrates the TC damagefunction with em-dat annual
%   damage data
%
% CALLING SEQUENCE:
%   result=calibrate_TC_DF_emdat_framework_region(TCBasinID,value_mode,cropped_assets,resolution,calibrate_countries,hazard_filename,number_free_parameters,years_considered)
%  
% EXAMPLE:
% TCbasins = {'CAR' 'NAM' 'NWP' 'NIN' 'SIN' 'PIS' 'AUS'};
% TCBasinIDs = [11    12    2     3     4     51    52];
% % ii =         1 .   2 .  3 .   4 .   5 .   6 .   7
% 
% ii_TCBasinID            = 5;%1:length(TCBasinIDs);
% value_mode              = 2;
% cropped_assets          = 0;
% resolution              = 300;
% calibrate_countries     = 0;
% hazard_filename         = 'GLB_0360as_TC_hist';
% number_free_parameters  = 2; 
% years_considered = 1980:2015;
% 
% for id = ii_TCBasinID   
%     disp(TCbasins{id});
%     
%     results_struct.(TCbasins{id})=calibrate_TC_DF_emdat_framework_region(...
%         TCBasinIDs(id),value_mode,cropped_assets,resolution,calibrate_countries,...
%         hazard_filename,number_free_parameters,years_considered);
% end
%
% EXAMPLE2:
% bsub -N -B -J "matlab_TC_cal003-2" -W 120:00 -R "rusage[mem=20000]" -oo logs/l_TC_cal_003-2.txt -e logs/e_TC_cal_003-2.txt 'matlab -nodisplay -nojvm -singleCompThread -r "calibrate_TC_DF_emdat_framework_region(2,2,0,300,0,'GLB_0360as_TC_hist_1000km.mat',2,1980:2015,1)"'
%
% INPUTS:
%   several (to be explained)
%
% OUTPUT:
%   result: resulting shape and scale parameters for damage function &
%   optimized function value
% 
% REQUIREMENTS:
%    
%      .../country_risk/data/NatID_RegID_basins_refined_201807.mat
%      .../country_risk/data/em_data file per region (e.g. em_data_TC_2005_CAR_1980-2015.mat)
%      .../entities/entity_region_file per region (for calibration)
%      required hazard file
% 
% MODIFICATION HISTORY:
% Samuel Eberenz, eberenz@posteo.eu, 20180718, init
% Samuel Eberenz, eberenz@posteo.eu, 20180905, add option to hand over entity file name
%-

if ~exist('on_cluster','var'), on_cluster=[];end
if isempty(on_cluster),on_cluster = 0;end
try, global climada_global; catch; end

if on_cluster, regions_from_xls = 0; end% xlsx can't be read on EULER;

if on_cluster && (~exist('climada_global','var') || isempty(climada_global))
    try
        cd ~/climada;
        startup      
    catch
        error('Location ~/climada not found. Please set location of climada folder on cluster correctly');
    end
end
clear hazard entity* x*
% https://ch.mathworks.com/help/gads/examples/constrained-minimization-using-pattern-search.html

%%

peril_ID = 'TC';

if ~exist('TCBasinID','var'),TCBasinID=[];end
if ~exist('value_mode','var'),value_mode=[];end
if ~exist('cropped_assets','var'),cropped_assets=[];end
if ~exist('resolution','var'),resolution=[];end
if ~exist('calibrate_countries','var'),calibrate_countries=[];end
if ~exist('hazard_filename','var'),hazard_filename=[];end
if ~exist('number_free_parameters','var'),number_free_parameters=[];end
if ~exist('years_considered','var'),years_considered=[];end
if ~exist('hand_over_entity_file','var'),hand_over_entity_file=[];end


% TCbasins = {'CAR' 'NAM' 'NWP' 'NIN' 'SIN' 'PIS' 'AUS'};
% TCBasinIDs = [11    12    2     3     4     51    52];
if isempty(TCBasinID),TCBasinID=1;end
if isempty(value_mode),value_mode=2;end
if isempty(cropped_assets),cropped_assets=0;end
if isempty(resolution),resolution=300;end
if isempty(calibrate_countries),calibrate_countries=0;end
if isempty(hazard_filename),hazard_filename=['GLB_0360as_',peril_ID,'_hist'];end
if isempty(number_free_parameters),number_free_parameters=2;end
if isempty(hand_over_entity_file),hand_over_entity_file=0;end


if ~exist('regions_from_xls','var'),regions_from_xls=0;end
if ~exist('reference_year','var'),reference_year=2005;end

delta_shape_parameter = 49; % v_half - v_threshold (m/s)
encode = 0;

optimizerType='R2';
remove_0_YDS_years = 1;
%optimizerType='R';
%optimizerType='logR';

full_parameter_search = 1;
fminconSwitch = 0;
save_output = 1;
force_overwrite_output = 0;

%

Input_path = [climada_global.modules_dir filesep 'country_risk' filesep 'data'];
savedir = [climada_global.results_dir filesep 'TC_calibration'];
if ~exist(savedir,'dir') && save_output
    system(['mkdir ' savedir]);
end

regions_file_xls = [Input_path filesep 'NatID_RegID_basins_refined_201807.xls'];
regions_file_mat = [Input_path filesep 'NatID_RegID_basins_refined_201807.mat'];

% get country / regions mapping
if regions_from_xls
    regions.countries = climada_xlsread(0,regions_file_xls,'NatID',1);
    regions.mapping = climada_xlsread(0,regions_file_xls,'TCBasins',1);
    save(regions_file_mat,'regions','-v7.3');
else
    load(regions_file_mat);
end
% Load EM-DAT data
reference_year_in = reference_year;
em_data_file_mat = [Input_path filesep 'em_data_' peril_ID '_' num2str(reference_year) '_' regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)} '_' num2str(years_considered(1)) '-' num2str(years_considered(end)) '.mat'];
load (em_data_file_mat) % includes em_data_region and country_list
if reference_year_in ~= reference_year, warning('em-data reference_year mismatch!');end
reference_year = reference_year_in;
clear reference_year_in

if isempty(years_considered),years_considered=em_data_region.year;
else % set em-data damage 0 for all years not member of years_considered:
    ix_years_considered=ismember(em_data_region.year,years_considered);
    em_data_region.damage_cal=em_data_region.damage_cal.*ix_years_considered;
    if calibrate_countries
        for i_c = 1:length(country_list)
        	em_data_region.(['damage_' country_list{i_c}])=em_data_region.(['damage_' country_list{i_c}]).*ix_years_considered;
        end
    end    
end
fprintf('A total of %i year(s) are considered for calibration (between %i and %i)\n',length(years_considered),min(years_considered),max(years_considered));
%%
switch value_mode
    case 1
        entity_region_file = ['TCBasin_' regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)} ...
            '_GDP_LitPop_BM2016_' num2str(resolution) 'arcsec_ry' num2str(reference_year)];
    case 2
        entity_region_file = ['TCBasin_' regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)}...
            '_NFW_LitPop_BM2016_' num2str(resolution) 'arcsec_ry' num2str(reference_year)];
    case 3      
        entity_region_file = ['TCBasin_' regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)}...
            '_TOW_LitPop_BM2016_' num2str(resolution) 'arcsec_ry' num2str(reference_year)];
end
if cropped_assets
	entity_region_file = [entity_region_file '_cropped'];
end
if exist([climada_global.entities_dir filesep entity_region_file '_cal.mat'],'file') % extra calibration entity?
    entity_region_file = [entity_region_file '_cal'];
end

% admin0_ISO3=country_list{i_country};
% [admin0_name,admin0_ISO3] = climada_country_name(admin0_ISO3); % get full name
% entity_filename = [admin0_ISO3 '_GDP_LitPop_BM2016_', num2str(resolution), 'arcsec_ry' num2str(reference_year)];
%
entity = climada_entity_load(entity_region_file);
hazard = climada_hazard_load(hazard_filename);
hazard = climada_hazard_reset_yearset(hazard,1);
if encode, entity = climada_assets_encode(entity,hazard,40000);end

if remove_0_YDS_years
    
    % figure; plot(em_data_region.year, em_data_region.damage_cal,'bo');
    
    v_threshold = 15;
    v_half = 20;
    scale = 1;
    choose_DF = strcmp(entity.damagefunctions.peril_ID,'TC') .* (entity.damagefunctions.DamageFunID==1);

    % set PPA to 1 for all intensities
    entity.damagefunctions.PAA(choose_DF==1) = 1;
    % set MDD with x(1)= v_half [m/s] % 
    % based on Elliott et al. 2015, similar also Sealy & Strobl et al. 2017 and Emanuel, 2011
    VV_temp = max((entity.damagefunctions.Intensity(choose_DF==1) - v_threshold),0)/(v_half-v_threshold);
    entity.damagefunctions.MDD(choose_DF==1)= VV_temp.^3 ./ (1+VV_temp.^3);
    % linearly damp MDD by multiplication with scale x(2), 0<x(2)<=1
    entity.damagefunctions.MDD(choose_DF==1)= scale * entity.damagefunctions.MDD(choose_DF==1);
    % (2) calculate YDS
    EDS = climada_EDS_calc(entity,hazard,[],[],2);
    [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
    [LIA,LOCB] = ismember(em_data_region.year,YDS.yyyy);
    % exclude years with 0-damage-entries:
    year_i_vector = find(LIA);
    year_i_vector(em_data_region.damage_cal(year_i_vector)==0)=[]; 
    year_i_vector(YDS.damage(LOCB(year_i_vector))==0)=[];
    em_data_region.damage_cal = em_data_region.damage_cal(year_i_vector);
    em_data_region.year = em_data_region.year(year_i_vector);
    for i_c = 1:length(country_list)
        em_data_region.(['damage_' country_list{i_c}])=em_data_region.(['damage_' country_list{i_c}])(year_i_vector);
    end
    % hold on;plot(em_data_region.year, em_data_region.damage_cal,'rx');hold off;
    N_years = length(find(em_data_region.damage_cal>0));
    disp(N_years)
end


%% set bounbdaries and starting values
% 
% TCbasins = {'CAR' 'NAM' 'NWP' 'NIN' 'SIN' 'PIS' 'AUS'};
% TCBasinIDs = [11    12    2     3     4     51    52];

% switch TCBasinID
%     case 11 % CAR / -
%         v_thres_0 = 25;
%     case 12 % NAM / USA
%         v_thres_0 = 29.8;
%     case 2 % NWP / JPN
%         v_thres_0 = 33.0;
%     case 3 % NIN / IND
%         v_thres_0 = 25.0;
%     case 4 % SIN / MDG
%         v_thres_0 = 25.0;
%     case 51 % PIS / -
%         v_thres_0 = 25.0;
%     case 52 % SWP / AUS
%         v_thres_0 = 25.0;
%    otherwise
        v_thres_0 = 30.7;
%end
        
        
switch number_free_parameters
    case 2
%         x0 = [.5*(max(v_thres_0,30) + v_thres_0) .5]; % starting values of free parameters
%         bounds.lb = [max(v_thres_0,30)-5. 1e-9]; % never set lower bound of scale to 0! Otherwise calibration goes wrong.
%         bounds.ub = [v_thres_0+5. 1.];        
        x0 = [v_thres_0 .5]; % starting values of free parameters
        bounds.lb = [max(v_thres_0-5,25.7) 1e-9]; % never set lower bound of scale to 0! Otherwise calibration goes wrong.
        bounds.ub = [v_thres_0+5. 1.];
    otherwise
        error('number_free_parameters other than 2 is not implemented yet');
end

norm.lb = ones(size(bounds.lb));
norm.ub = norm.lb+1;

norm.x0 = (x0-bounds.lb)./(bounds.ub-bounds.lb) .* (norm.ub-norm.lb) + norm.lb;


%%
% figure; climada_entity_plot(entity,3);
if hand_over_entity_file
    clear entity;
    entity = entity_region_file;
end

% define anonymous function with input factor x (parameters of the damage
% function):
fun = @(x)calibrate_TC_DF_emdat_region(x,delta_shape_parameter, entity, hazard, em_data_region, norm, bounds,optimizerType); % sets all inputvar for the function except for x, use normalized x.

if fminconSwitch
%     There are no linear constraints, so set those arguments to []:
%     [x,fval,exitflag,output]=fmincon(fun,norm.x0,[],[],[],[],norm.lb,norm.ub,[],fmincomOptions)
%     Try an initial point in the middle of the region. Find the minimum of fun, subject to the bound constraints.
    fmincomOptions = optimoptions('fmincon','MaxFunctionEvaluations',80,'Display','iter','StepTolerance',0.001);
    results.x=[];
    results.fval=[];
    results.exitflag=[];
    results.ME=[];
    
    try
        tic
        [results.x,results.fval,results.exitflag]=...
            fmincon(fun,norm.x0,[],[],[],[],norm.lb,norm.ub,[],fmincomOptions);
        toc
        results.x=(results.x-norm.lb).*(bounds.ub-bounds.lb)./(norm.ub-norm.lb)+bounds.lb;
        results.fval
    catch ME
        warning(['fmincon failed. Check results.ME for more info.']);
        display(ME.identifier)
        display(ME.message)
        results.ME = ME;
        clear ME
    end

    results.resolution = resolution;
    if save_output
        save_file_name=[savedir filesep regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)}...
            '_' peril_ID '_calibrate_region_fmincon_litpop_gdp_' num2str(number_free_parameters) '-' num2str(resolution) '.mat'];
        
        while (~force_overwrite_output && exist(save_file_name,'file')),save_file_name=strrep(save_file_name,'.mat','_.mat');end % avoid overwriting
        
        save(save_file_name,'TCBasinID','peril_ID','x0','results','years_considered','cropped_assets','resolution','value_mode','-v7.3');
    end
end

%parpool('local_small')
if full_parameter_search
    
    options = optimoptions('patternsearch','UseParallel',false,...
    'UseCompletePoll', true, 'UseVectorized', false,...
    'MaxFunctionEvaluations',1200,'Display','iter',...
    'Cache','on','InitialMeshSize',.25,...
    'PollMethod','GPSPositiveBasis2N','StepTolerance',0.001);
    tic
    [x_result,fval] = patternsearch(fun,norm.x0,[],[],[],[],norm.lb,norm.ub,[],options);
    toc
    result.region=(x_result-norm.lb).*(bounds.ub-bounds.lb)./(norm.ub-norm.lb)+bounds.lb
    fval
    if save_output && ~calibrate_countries
        
        save_file_name=[savedir filesep regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)}...
            '_' peril_ID '_decay_region_calibrate_litpop_gdp_' num2str(number_free_parameters) '-' num2str(value_mode) '-' num2str(resolution) '.mat'];
        
        while (~force_overwrite_output && exist(save_file_name,'file')),save_file_name=strrep(save_file_name,'.mat',['_' hazard_filename '_' datestr(today) '.mat']);end % avoid overwriting
        
        save(save_file_name,'result','fval','hazard_filename','TCBasinID');
        
    end
end
% clear entity hazard

% calibration of single countries with more than or equal to N_min damage years
if calibrate_countries
    N_min = 1;
    
    clear entity
    em_data_country.year = em_data_region.year;
    result.result_c = NaN*ones(length(country_list),length(x0));
    result.fval_c = NaN*ones(length(country_list),1);
    N_damageyears = result.fval_c;
    options = optimoptions('patternsearch','UseParallel',false,...
        'UseCompletePoll', true, 'UseVectorized', false,...
        'MaxFunctionEvaluations',1200,'Display','iter',...
        'Cache','on','InitialMeshSize',.25,...
        'PollMethod','GPSPositiveBasis2N','StepTolerance',0.001);
    
    for i_c = 1:length(country_list)
        %%
        N_damageyears(i_c) = sum(em_data_region.(['damage_' country_list{i_c}])>0);
        disp([country_list{i_c} ', N=' num2str(N_damageyears(i_c))]);
        if N_damageyears(i_c) >= N_min
            entity = climada_entity_load([country_list{i_c} '_GDP_LitPop_BM2016_' num2str(resolution) 'arcsec_ry' num2str(reference_year) '.mat']);
            em_data_country.damage = em_data_region.(['damage_' country_list{i_c}]);
            fun = @(x)calibrate_TC_DF_emdat_region(x,delta_shape_parameter, entity, hazard, em_data_country, norm, bounds,optimizerType);
            if sum(entity.assets.Value)>0
                [x_result,result.fval_c(i_c)] = patternsearch(fun,norm.x0,[],[],[],[],norm.lb,norm.ub,[],options);
            
                result.result_c(i_c,:)=(x_result-norm.lb).*(bounds.ub-bounds.lb)./(norm.ub-norm.lb)+bounds.lb
            else
                %result.result_c(i_c,:)=[NaN NaN];
                warning('entity contains no values');
            end
        end
    end
    %%
    if save_output
        save_file_name = [savedir filesep regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)}...
            '_' peril_ID '_decay_region_countries_calibrate_litpop_gdp_' num2str(number_free_parameters) '-' num2str(value_mode) '-' num2str(resolution) '.mat'];
        
        while (~force_overwrite_output && exist(save_file_name,'file')),save_file_name=strrep(save_file_name,'.mat','_.mat');end % avoid overwriting
        
        save(save_file_name,'country_list','result','fval','resolution','years_considered','-v7.3');
        save_file_name
    end
end


if on_cluster, cd /cluster/home/eberenzs/ ;exit; end;


end
