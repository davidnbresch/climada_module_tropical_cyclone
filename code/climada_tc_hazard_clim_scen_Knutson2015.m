function hazard = climada_tc_hazard_clim_scen_Knutson2015(hazard,tc_tracks,target_rcp_scenario,target_year,make_plots,output_filename)
% MODULE:
% tropical_cyclone
% NAME:
% climada_tc_hazard_clim_scen_Knutson2015
% PURPOSE:
%   For different years and RCPs, compute the future TC hazard set based on: 
%   'Knutson2015' :: Knutson et al., 2015. Global projections of intense tropical 
%                    cyclone activity for the late twenty-first century from dynamical
%                    downscaling of CMIP5/RCP4.5 scenarios.
%                    (Comparison 2081-2100 (i.e., late twenty-first
%                    century) and 2001-20 (i.e., present day))
%                    Late twenty-first century effects on intensity and
%                    frequency per Saffir-Simpson-category and ocean basin
%                    is scaled to target year and target RCP proportional
%                    to total radiative forcing of the respective RCP and
%                    year.
% CALLING SEQUENCE:
%   hazard = load(hazard_file); % OR hazard='hazard_file_name';
%   tc_tracks = tc_track=climada_tc_track_load(tc_track_filename,0);
%   hazard = climada_tc_hazard_clim_scen_Knutson2015(hazard,tc_tracks,target_rcp_scenario,target_year,make_plots);
% EXAMPLES:
%   hazard = climada_tc_hazard_clim_scen_Knutson2015([],[],45,2045,1,'AUTO'); % 2-degree, 2045
%   hazard = climada_tc_hazard_clim_scen_Knutson2015([],[],85,2045,1,'AUTO'); % 4-degree, 2045
%
%   hazard = climada_tc_hazard_clim_scen_Knutson2015(hazard,tc_tracks,60,2050,1,'NO_SAVE');
%   
%   hazard = climada_tc_hazard_clim_scen_Knutson2015([],[],45,2080,0,'AUTO');
%
%   hazard = climada_tc_hazard_clim_scen_Knutson2015('GLB_0360as_TC_hist',[climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs.mat'],60,2030,1);
%   
% OPTIONAL INPUT PARAMETERS:
%   hazard:     hazard set today (TC) as file name or struct. default: 'GLB_0360as_TC'
%   centroids:  cetroids of hazard set as file or struct (only required if
%       hazard.basins is missing) default: [climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs_prob.mat']
%   target_rcp_scenario: number of CMIP5/RCP scenario without decimal point, i.e. 26 for
%       RCP2.6. Possible Inputs: 26, 45 (default), 60, 85.
%   target_year: target year of future hazard. default=2050
%   make_plots: make plots? default=0
%   output_filename: the name of the newly TC hazard event set
%      (if ='AUTO' the hazard file name is derived automatically)
%      (if ='NO_SAVE' or ='NOSAVE', the hazard is just returned, not saved)
%       > default 'NO_SAVE' if not given.
% OUTPUTS:
%   hazard:     new TC hazard for a future reference year, given a climate
%               change scenario specified by screw.
%
% MODIFICATION HISTORY:
% Samuel Eberenz, eberenz@posteo.eu, 20180705, init
% david.bresch@gmail.com, 20180815, PARAMETERS label introduced and parameters in code grouped, output filename streamlined
% Samuel Eberenz, eberenz@posteo.eu, 20180705, debug ylim and xlim if change in f or intensity is 0
% david.bresch@gmail.com, 20190426, output_filename fixed
%-
%% initiate

global climada_global;
if ~climada_init_vars,return;end % init/import global variables

if ~exist('hazard'    , 'var'), hazard     = []; end
if ~exist('tc_tracks' , 'var'), tc_tracks  = []; end
if ~exist('target_rcp_scenario' , 'var'), target_rcp_scenario  = []; end
if ~exist('target_year' , 'var'), target_year  = []; end
if ~exist('make_plots'  , 'var'), make_plots   = []; end
if ~exist('output_filename', 'var'), output_filename = []; end

% PARAMETERS
%
if isempty(hazard)
    hazard='GLB_0360as_TC';
    fprintf('starting from default %s hazard set\n',hazard);
end
if isempty(tc_tracks), tc_tracks=[climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs_prob.mat'];end
if isempty(target_rcp_scenario), target_rcp_scenario=45; end
if isempty(target_year), target_year=2050; end
if isempty(make_plots), make_plots= 0; end
if isempty(output_filename), output_filename = 'NO_SAVE'; end
%
if contains(output_filename,'NOSAVE') || contains(output_filename,'NO_SAVE')
    save_output = 0;
elseif contains(output_filename,'AUTO')
    save_output = 1;
    %output_filename = sprintf('TC_hazard_Knutson2015_rcp%i_%i',target_rcp_scenario,target_year);
    output_filename = sprintf('TC_rcp%i_%i',target_rcp_scenario,target_year); % default generic name
    if ischar(hazard)
        hazard_filename_front = strsplit(hazard,'.mat'); % cell
        %output_filename = sprintf('%s_Knutson2015_rcp%i_%i',hazard_filename_front{1},target_rcp_scenario,target_year);
        output_filename = sprintf('%s_rcp%i_%i',hazard_filename_front{1},target_rcp_scenario,target_year);
    elseif isstruct(hazard)
        if isfield(hazard,'filename')
            [fP,fN,fE]=fileparts(hazard.filename);
            output_filename = sprintf('%s%s%s_rcp%i_%i%s',fP,filesep,fN,target_rcp_scenario,target_year,fE);
        end
    end
else
    save_output = 1;
end
%
save4octave = 0;
if climada_global.octave_mode,save4octave=1;end
%
% TC_TS = 0;
% elevation_data = 0; 
% centroid_filename=[climada_global.centroids_dir filesep 'GLB_NatID_grid_0360as_adv_1.mat'];
%
% chose reference (currently only Knutson 2015 implemented!)
knutson2015=1;
%
basins = {'NA' 'EP' 'WP' 'NI' 'SI' 'SP'}; %% pick basins.... 
%                         'GL' :: Global (default)
%                         'EP' :: North East Pacific Ocean
%                         'NA' :: North Atlantic Ocean
%                         'NI' :: North Indian Ocean
%                         'SA' :: South Atlantic Ocean
%                         'SI' :: South Indian Ocean
%                         'SP' :: South Pacific Ocean
%                         'WP' :: North West Pacific Ocean
filename_extension = '_all_basins_';


if knutson2015
    screw_ref = 'Knutson2015';
    % get time scale:
    time_scale=climada_screw_set_time_scaling(target_rcp_scenario,target_year,make_plots,45,[2001 2020],[2081 2100]);
else
    screw_ref = 'IPCC_SREX';
    warning('Not yet implemented properly for IPCC_SREX screws, time scale might be off.')
    time_scale=climada_screw_set_time_scaling(target_rcp_scenario,target_year,make_plots,45,[2001 2020],[2081 2100]);
end
disp(['time_scale = ' num2str(time_scale)]);
if time_scale>1
    warning('Radiative forcing of selected scenario is higher than forcing in RCP4.5 end of 21st century (reference scenario). Changes in intensity & frequency are scaled with a factor larger than 1.')
end
% check whether hazard is filename or struct, get hazard and basins:
if  ischar(hazard) || isstr(hazard)
    hazard_filename_front = strsplit(output_filename,'.mat'); % cell
    hazard=climada_hazard_load(hazard); 
else
    hazard_filename_front = ['TC_hazard_' datestr(now)];
end
if ~isfield(hazard,'basin') || (length(hazard.basin) ~= length(hazard.event_ID))
    if ischar(tc_tracks) || isstr(tc_tracks)  
        tc_tracks=climada_tc_track_load(tc_tracks,0);
    end
    hazard.basin = extractfield(tc_tracks,'basin');
    clear tc_tracks
end
% check consistency of hazard and tc_tracks
if length(hazard.basin) ~= length(hazard.event_ID), error('basin data from tc_tracks do not match hazard set');end

% change intensity and frequency of hazard:
screw = climada_hazard_climate_set_screw(screw_ref,'TC',basins);

hazard_cc = climada_hazard_climate_screw(hazard, 'NO_SAVE', target_year, screw, time_scale);

% set negative values 0:
hazard_cc.frequency(hazard_cc.frequency<0)=0;
hazard_cc.intensity(hazard_cc.intensity<0)=0;
% add comment to hazard file:
switch target_rcp_scenario
    case 26
        path_str = 'RCP2.6';
    case 45
        path_str = '2-deg path (RCP4.5)';
    case 60
        path_str = 'RCP6.0';
    case 85
        path_str = '4-deg path (RCP8.5)';
    otherwise
        path_str ='';
end
routine = 'climada_tc_hazard_clim_scen_Knutson2015';
hazard_cc.comment_clim=sprintf('Climate Change scenario %i, %s, adding basin-specific increments, see %s, %s',...
    target_year,path_str,routine,datestr(now));
disp(hazard_cc.comment_clim);

if time_scale>1
    hazard_cc.warning = 'Radiative forcing of selected scenario is higher than forcing in RCP4.5 end of 21st century (reference scenario). Changes in intensity & frequency are scaled with a factor larger than 1.';
end

% %% TC --> TS: Applying climada_ts_hazard_set.m for global grid TC hazard
% if TC_TS % not fully tested, not included in options yet
%     
%     close all;
%     if TC_screw, TC_TS, hazard4TS_filename = output_filename;
%     else, hazard4TS_filename = hazard_filename_front; end
% 
%     hazard_TS_filename=replace(hazard4TS_filename,'TC','TS');
% 
%     hazard=climada_hazard_load(hazard4TS_filename);  
%     if ~save_output
%         hazard_TS_filename='NO_SAVE';
%     end
%     hazard_cc=climada_ts_hazard_set(hazard,hazard_TS_filename,elevation_data,make_plots);
%     
%         
%     if save_output && save4octave, climada_save_mat_for_octave(hazard_TS_filename,'hazards',6,1,1);end
% end

if make_plots
    
    if ~exist('hazard_cc','var')
        hazard_cc = climada_hazard_load(hazard_cc_file);
    end
    if ~exist('hazard','var')
        hazard=climada_hazard_load(hazard_filename);
    end
    hazard_delta = hazard_cc;
    
    hazard_delta.intensity = hazard_cc.intensity - hazard.intensity;
    hazard_delta.frequency = hazard_cc.frequency - hazard.frequency;
    randomi= randperm(size(hazard_delta.intensity,2),1000);
    %%
    basins_all = unique(hazard_cc.basin);
    for i = 1:length(basins_all)
        IndexC = strfind(hazard_cc.basin, basins_all{i});
        Index = find(not(cellfun('isempty', IndexC)));
        diff_freq(i) = nanmean(nanmean(hazard_delta.frequency(Index)));
        INTENSITY = hazard_delta.intensity(Index,randomi)./hazard.intensity(Index,randomi);
        for j=1:size(INTENSITY,1)
            DIFF_INT(j) = (nanmean(INTENSITY(j,:)));
        end
        diff_int(i) = nanmean(DIFF_INT(:));
        clear Index* INTENSITY DIFF_INT
    end

    figure; 
    yyaxis left
    ylabel('Frequency')
    stem(diff_freq./hazard.frequency(1),'*','MarkerSize',10);
    max_change = max(abs(diff_freq./hazard.frequency(1)));
    if max_change == 0, max_change = 1;end
    ylim([-1.1*max_change 1.1*max_change])
    yyaxis right
    ylabel('Intensity')
    stem(diff_int,'d','MarkerSize',10)
    max_change = max(abs(diff_int));
    if max_change == 0, max_change = 1;end
    ylim([-1.1*max_change 1.1*max_change]);
    xticklabels(basins_all);
    title('Tropical Cyclones: Relative Delta per Basin');
    xlabel('Ocean Basin')
    hold on
    plot(0:8,zeros(9,1),':k');
    hold off
    grid on
    legend('Frequency','Intensity')
    params.difference=1;
    figure;res = climada_hazard_plot(hazard_delta,0,[],params);set(gcf, 'Position', [100, 100, 1100, 600])
    title('Difference in maximum intensity between output and input hazard')
%    figure;res = climada_hazard_plot(hazard_cc,0);set(gcf, 'Position', [300, 300, 1100, 600])

%     
%     i = 2;%1:length(basins_all)
%     IndexC = strfind(hazard_cc.basin, basins_all{i});
%     Index_only = find(not(cellfun('isempty', IndexC)));
%     Index_exclude = find((cellfun('isempty', IndexC)));
%     hazard_exclude = hazard_cc;
%     hazard_only = hazard_cc;
%     hazard_exclude.intensity(Index_only,:)=0;
%     hazard_only.intensity(Index_exclude,:)=0;
%     figure; climada_hazard_plot(hazard_only,0); title('only')
%     figure; climada_hazard_plot(hazard_exclude,0); title('exclude')
%     clear Index* INTENSITY DIFF_INT
    
end

hazard = hazard_cc;
if save_output 
    save(output_filename,'hazard')
    cprintf([113 198 113]/255, 'climate change scenario hazard set saved in %s\n',output_filename);
    if save4octave, climada_save_mat_for_octave(output_filename,'hazards',6,1,1);end
end



end




