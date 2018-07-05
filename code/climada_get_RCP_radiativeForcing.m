function [rcp_RF, rcp_years] = climada_get_RCP_radiativeForcing(rcp_scenario,years,make_plot,rcp_file)
%%%
%   
% MODULE:
%   tropical_cyclone
% NAME:
%   climada_get_RCP_radiativeForcing
% PURPOSE:
%   Return total radiative forcing for selected RCP
%   Note: Requires rcp_db.xls from http://www.iiasa.ac.at/web-apps/tnt/RcpDb > Compare > Climate indicators > Radiative forcing > Total  

% EXAMPLES:
%   [rcp85_RF, years] = climada_get_RCP_radiativeForcing(85,[],1);
%   
%   climada_get_RCP_radiativeForcing(26,2050,0);  
%
% INPUT (OPTIONAL):
% rcp_scenario: number of RCP scenario without decimal point, i.e. 26 for
%       RCP2.6. Possible Inputs: 26, 45 (default), 60, 85.
% years: integer or vector with year(s) for which results are to be returned (default 2000:5:2100)
% make_plot: boolean, make a plot of all four RCPs (default=0)
% rcp_file: user defined path of RCP file (XLS), default = .../climada_modules/tropical_cyclone/data/rcp_db.xls
%
% OUTPUT:
%   rcp_RF: vector with value(s) of total radiative forcing in W/m^2 for chosen RCP and years
%   rcp_years: integer or vector with years

% MODIFICATION HISTORY:
%  Samuel Eberenz, eberenz@posteo.eu, 20180704, initial.
%

%%
% init variables

global climada_global;
if ~climada_init_vars,return;end % init/import global variables

%%

rcp_RF = [];
if ~exist('rcp_scenario','var'),rcp_scenario=[];end
if ~exist('years','var'),years=[];end
if ~exist('make_plot','var'),make_plot=[];end
if ~exist('rcp_file','var'),rcp_file=[];end

if isempty(rcp_scenario),rcp_scenario=45;end
if isempty(make_plot),make_plot=0;end
if isempty(rcp_file),rcp_file=[climada_global.modules_dir filesep 'tropical_cyclone' filesep 'data' filesep 'rcp_db.xls'];end
%%
if exist(rcp_file,'file')
    rcp_RF = xlsread(rcp_file);
    rcp_RF_import_headers = climada_xlsread(0,rcp_file,[],1);
    rcp_years = rcp_RF(1,5:end);
    rcp_RF = rcp_RF(2:end,5:end);
    rcp_scenarios = [60 45 26 85]; % make sure order of scenarios is the same as in rcp_file
        
    if make_plot
        colors={'m','b','g','r','k','y'};
        
        figure();
        set(gca,'FontSize',15)%,'fontWeight','bold')
        hold on
        for i=1:length(rcp_scenarios)
            plot(rcp_years,rcp_RF(i,:),colors{i},'LineWidth',3);
        end
        legend(rcp_RF_import_headers.Scenario{1:i},'Location','NorthWest')
        xlabel('Year')
        ylabel('Total Radiative Forcing [W/m^2]')
        title('Radiative Forcing CMIP5 RCPs')
        
        grid on
    end
    
    i_scenario = find(rcp_scenario==rcp_scenarios);
    rcp_RF = rcp_RF(i_scenario,:);
    if isempty(i_scenario)
        warning('RCP Scenario misspecified. Try 26 45 60 or 85.')
    end
    
    if ~isempty(years) % interpolate if years are given
        if min(years)<2000 || max(years)>2100
            warning('Requested years outside period 2000 to 2100')
        end
        
        rcp_RF = interp1(rcp_years,rcp_RF,years);
        rcp_years=years;
        if make_plot
            plot(years,rcp_RF,'xk','MarkerSize',15,'LineWidth',3)
            hold off
            legend(rcp_RF_import_headers.Scenario{1:i},'Selected Forcing','Location','NorthWest')
        end
        
    end
    
else
    warning(['RCP File (' rcp_file ') not found. Download from http://www.iiasa.ac.at/web-apps/tnt/RcpDb > Compare > Climate indicators > Radiative forcing > Total '])
    return
end

end
