function time_scale = climada_screw_set_time_scaling(rcp_scenario,target_years,make_plot,screw_rcp_scenario,screw_base_years,screw_end_years)
% MODULE:
%   tropical_cyclone
% NAME:
%   climada_screw_set_time_scaling
% PURPOSE:
%   get a scaling factor for climate change screws (i.e. for TC) which is 
%   essentially the ratio between the total radiative forcing of a chosen
%   RCP-scenario for (a) certain year(s) and the radiative forcing and time
%   period assumed in a particular "screw"
%
% CALLING SEQUENCE:
%  screw = climada_hazard_climate_set_screw(reference, hazard_type, basin)
%  time_scale = climada_screw_set_time_scaling(45,reference_year)
%  hazard = climada_hazard_climate_screw(hazard, hazard_set_file, reference_year, screw, time_scale)
%
% EXAMPLES:
%   time_scale = climada_screw_set_time_scaling(60,2050,1,60,2000,2100)
%   % or:
%   target_year = 2050;
%   screw = climada_hazard_climate_set_screw('Knutson2015, 'TC', 'WP')
%   time_scale = climada_screw_set_time_scaling(85,target_year,0,45,2001:2020,2081:2100)
%   hazard = climada_hazard_climate_screw(hazard, hazard_set_file, target_year, screw, time_scale)
%
% INPUTS:
%   rcp_scenario: number of RCP scenario without decimal point, i.e. 26 for
%       RCP2.6. Possible Inputs: 26, 45, 60, 85.
%   target_years: integer or vector with year(s) for which scaling factor is to be returned

% OPTIONAL:
%   make_plot: boolean, make a plot of all four RCPs' radiative forcings (default=0)
%   screw_rcp_scenario:  number of RCP scenario OF SCREW without decimal point, i.e. 26 for
%       RCP2.6. Possible Inputs: 26, 45 (default), 60, 85.
%   screw_base_years: vector of years that are the baseline period (today's
%       climate) in the screw, default: [2001 2020]
%   screw_end_years: vector of years that are the target period (future's
%       climate) in the screw, default: [2081 2100]
% OUTPUTS:
%   time_scale: double with a scaling factor for the screw based on
%       RCP radiative forcing ratios
% MODIFICATION HISTORY:
%  Samuel Eberenz, eberenz@posteo.eu, 20180705, initial.    
%%%

% initiate
time_scale = [];
if ~exist('rcp_scenario','var'),return;end % 26 45 60 or 85
if ~exist('target_years','var'),return;end

if ~exist('make_plot','var'),make_plot=0;end % Knutson 2015 defaults
if ~exist('screw_rcp_scenario','var'),screw_rcp_scenario=45;end % Knutson 2015 defaults
if ~exist('screw_base_years','var'),screw_base_years=[2001 2020];end % Knutson 2015 defaults
if ~exist('screw_end_years','var'),screw_end_years=[2081 2100];end % Knutson 2015 defaults

% Get radiative forcings (RF): climada_get_RCP_radiativeForcing(rcp_scenario,years,make_plot,rcp_file)
RF_base = nanmean(climada_get_RCP_radiativeForcing(screw_rcp_scenario,screw_base_years(1):1:screw_base_years(end)));
RF_end  = nanmean(climada_get_RCP_radiativeForcing(screw_rcp_scenario,screw_end_years(1):1:screw_end_years(end)));
RF_target = nanmean(climada_get_RCP_radiativeForcing(rcp_scenario,target_years(1):1:target_years(end),make_plot));

% compute scale as ratio of forcings
time_scale = max((RF_target-RF_base)/(RF_end-RF_base),0);
    
end
