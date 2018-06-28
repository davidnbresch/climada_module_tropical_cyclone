function damagefunctions=climada_tc_damagefun_emanuel2011(damagefunctions,v_threshold,v_half,scale,DamageFunID,intensity_step_size,max_intensity,plot_df)
% MODULE:
%   tropical_cyclome
% NAME:
%   climada_inshape_new
% PURPOSE:
%   Create a TC damage function (mapping MDD on wind speed intensity) based
%   on the cubic damage function shape proposed by: 
%       Emanuel, Kerry. ?Global Warming Effects on U.S. Hurricane Damage.? 
%       Weather, Climate, and Society 3, no. 4 (October 2011): 261?68. 
%       https://doi.org/10.1175/WCAS-D-11-00007.1.
%   The damage function is written to a new or provided damagefunction-struct,
%   i.e. an entity.damagefunctions.
%   PAA is set to 1 in this function.
%
% CALLING SEQUENCE:
% EXAMPLES:
%   entity.damagefunctions = climada_tc_damagefun_emanuel2011(entity.damagefunctions,25,60,0.2);
%
%   new_damagefunction1 = climada_tc_damagefun_emanuel2011([],20,50,0.5);
%
%   new_damagefunction2 = climada_tc_damagefun_emanuel2011();
%
% INPUTS:
% ALL OPTIONAL:
%   damagefunctions: existing damagefunctions struct to append with new damagefunction
%   v_threshold: First shape parameter, wind speed in m/s below which no
%       damage is expected. Default = 25 m/s
%   v_half: Second shape parameter,wind speed in m/s at which 50% of max.
%       damage is expected. Default = v_threshold + 49 m/s.
%   scale: Scale parameter, linear scaling of MDD. 0<=scale<=1. Default = 1. 
%   DamageFunID: int, ID of new damagefunction. Default = 1.
%   intensity_step_size: resolution of intensity in m/s. Default = 5 m/s
%   max_intensity: maximum intensity in m/s. Default = 120 m/s
%   plot_df: boolean, plot new damagefunction? default = 0
% OUTPUTS:
%   damagefunctions: struct with new damagefunction or appended damagefunctions.
% MODIFICATION HISTORY:
%  Samuel Eberenz, eberenz@posteo.eu, 20180628, initial.
%
%%

global climada_global;
if ~climada_init_vars,return;end % init/import global variables
delta_shape_parameter_default = 49; % m/s, for calculation of v_half if not provided,

if ~exist('damagefunctions','var'),damagefunctions=[];end
if ~exist('scale','var'),scale=[];end
if ~exist('v_threshold','var'),v_threshold=[];end
if ~exist('v_half','var'),v_half=[];end
if ~exist('DamageFunID','var'),DamageFunID=[];end
if ~exist('intensity_step_size','var'),intensity_step_size=[];end
if ~exist('max_intensity','var'),max_intensity=[];end
if ~exist('plot_df','var'),plot_df=[];end

if isempty(scale),scale=1;end % m/s
if isempty(v_threshold),v_threshold=25;end % m/s
if isempty(v_half) 
    v_half = v_threshold+delta_shape_parameter_default;
    fprintf('v_half not provided; v_half was set to %0.1f m/s + %0.1f m/s = %0.1f m/s.\n',v_threshold,delta_shape_parameter_default,v_half);
end 

if isempty(intensity_step_size),intensity_step_size=5;end % m/s
if isempty(max_intensity),max_intensity=120;end % m/s
if isempty(DamageFunID),DamageFunID=1;end 
if isempty(plot_df),plot_df=0;end 


if v_half <= v_threshold || v_threshold <0 || scale <=0
    error('shape or scale parameter out of range')
end
if scale > 1
    warning('scale parameter larger than 1.0')
end

name = sprintf('Tropical Cyclone Emanuel %0.0f %0.0f %0.2f', v_threshold, v_half, scale);

% initiate new damagefunction struct:
damagefunction_new.filename='';
damagefunction_new.Intensity=0:intensity_step_size:max_intensity;
damagefunction_new.DamageFunID = zeros(size(damagefunction_new.Intensity))+DamageFunID;
damagefunction_new.PAA = ones(size(damagefunction_new.Intensity)); % PPA set to 1 for all intensities
damagefunction_new.MDD = zeros(size(damagefunction_new.Intensity));
damagefunction_new.peril_ID = repmat({'TC'},size(damagefunction_new.Intensity));
damagefunction_new.Intensity_unit = repmat({'m/s'},size(damagefunction_new.Intensity));
damagefunction_new.name = repmat({name},size(damagefunction_new.Intensity));
damagefunction_new.datenum = zeros(size(damagefunction_new.Intensity))+ now;

% Set MDD according to shape and scale parameters:

v_temporary = max((damagefunction_new.Intensity-v_threshold),0)/(v_half-v_threshold);
damagefunction_new.MDD= v_temporary.^3 ./ (1+v_temporary.^3);
% linearly dampen MDD by multiplication with scale x(2), 0<x(2)<=1
damagefunction_new.MDD= scale * damagefunction_new.MDD;
clear v_temporary 
if plot_df;figure;climada_damagefunctions_plot(damagefunction_new,'TC');end

%% append damagefunctions with new damagefunction
if isempty(damagefunctions)
    damagefunctions = damagefunction_new;
else
    choose_DF = strcmp(damagefunctions.peril_ID,'TC') .* (damagefunctions.DamageFunID==DamageFunID);
    if max(choose_DF)
        
        fprintf('WARNING: original TC %03.0f was reassigned to TC %03.0f.\n', DamageFunID, max(damagefunctions.DamageFunID)+1);
        damagefunctions.DamageFunID(find(choose_DF))=max(damagefunctions.DamageFunID)+1;
    end
    %%
    damagefunctions_merged = struct;  %final structure
    for field = fieldnames(damagefunctions)'
        fname = field{1};
        damagefunctions_merged.(fname) = [damagefunction_new.(fname) damagefunctions.(fname)];
    end
    damagefunctions=damagefunctions_merged;
    clear damagefunctions_merged damagefunction_new
end

end
