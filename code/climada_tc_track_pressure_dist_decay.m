function tc_track = climada_tc_track_pressure_dist_decay(tc_track, p_rel, check_plot)

% NAME:
%   climada_tc_track_pressure_dist_decay
% PURPOSE:
%   Incorporate pressure decay after landfall for probabilistic tracks, based on
%   historical tc tracks. The
%   decay can be calculated based on historical tracks or can be loaded
%   from a mat file (p_rel)
%   for:  isimip_windfield_holland
% CALLING SEQUENCE:
%     [~,p_rel_v]  = climada_tc_track_wind_decay_calculate(tc_track,check_plots);
%     tc_track_decay = climada_tc_track_wind_decay(tc_track_prob, p_rel_v, check_plots);
%     [~,p_rel_p]  = climada_tc_track_pressure_decay_calculate(tc_track,check_plots);
%     tc_track_decay = climada_tc_track_pressure_decay(tc_track_decay, p_rel_p, check_plots);
%     hazard   = isimip_tc_hazard_set(tc_track_decay,hazard_file,centroids);
%     hazard = climada_hazard_reset_yearset(hazard,1);
%     save([climada_global.hazards_dir filesep 'hazard_file'],'hazard','-v7.3');
% EXAMPLE:
%   tc_track = climada_tc_track_pressure_decay(tc_track)
% INPUTS:
%   none, if tc_track empty prompted for
% OPTIONAL INPUT PARAMETERS:
%   p_rel: parameters for exponential decay, 
%          pressure decay = S-(S-1)*exp(-A*x), where A = p_rel(1,1), S=EnvironmentalPressure/PressureAtLandfall=p_rel(1,2)
%          A = p_rel(1,1), S = p_rel(1,2), can be newly calculated or can
%          be loaded from data within globalGDP modul
%   check_plot: to create plot
% OUTPUTS:
%   same structure now CentralPressure in tc track decays after landfall
% RESTRICTIONS:
%   changes CentralPressure only.
% MODIFICATION HISTORY:
% Samuel Eberenz, eberenz@posteo.eu, 20180713, adapted from climada_tc_track_wind_decay
%               to change CentralPressure p instead of MaxSustainedWind v
%               and relative to distance on land instead of time after landfall
% Samuel Eberenz, eberenz@posteo.eu, 20180716,changed definition of S, correct for CentralPressure > EnvironmentalPressure
% Samuel Eberenz, eberenz@posteo.eu, 20180717,reset definition of S = tc_track(t_i).EnvironmentalPressure(end)/min(tc_track(t_i).EnvironmentalPressure(end),p_landfall);
%-

% init global variables
global climada_global
if ~climada_init_vars, return; end

% check inputs, and set default values
if ~exist('tc_track'       , 'var'), tc_track      = []  ; end
if ~exist('p_rel'          , 'var'), p_rel         = []  ; end
if ~exist('check_plot'     , 'var'), check_plot    = 1   ; end

% prompt for tc_track if not given
if isempty(tc_track)
    tc_track             = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    tc_track_default     = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select tc track .mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select tc tracks:',tc_track_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tc_track = fullfile(pathname,filename);
    end
end
% load the tc track set, if a filename has been passed
if ~isstruct(tc_track)
    tc_track_file = tc_track;
    tc_track      = [];
    load(tc_track_file);
end
%%
% refine tc_tracks to 1 h
tc_track = climada_tc_equal_timestep(tc_track,1);

% make a copy
tc_track_ori_1h = tc_track;

if ~isfield(tc_track,'onLand')
    % find nodes on land and over sea
    tc_track = climada_tc_on_land(tc_track);
end
if ~isfield(tc_track,'distOnLand_km')
    % find nodes on land and over sea, and distance over land (with distance dependent recovery over water)
    tc_track = climada_tc_track_distanceOnLand(tc_track); 
end
%%

% number of generated and historical tracks
no_ori = sum([tc_track(:).orig_event_flag]);
no_gen = length(tc_track)/no_ori;

if isempty(p_rel)
    % try to load p_rel as variable from data folder
    modul_data_dir = [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
    filename       = [modul_data_dir filesep 'wind_decay_after_landfall_p_rel.mat'];
    A = exist(filename,'file');
    if A==2
        load(filename)
    end
end

v_scale_kn   = [34 64 83 96 113 135 500];
no_cat       = size(v_scale_kn,2);
p_vs_lf_time = cell(1,no_cat);
dist_vs_lf_time = p_vs_lf_time;
cmap         = jet(no_cat);

%% calculate exponential decay of wind speed after landfall if not given
if isempty(p_rel)
    fprintf('Calculate p rel (parameters for exponential wind decay) based on historical tracks...\n')
    [~, p_rel] = climada_tc_track_pressure_dist_decay_calculate(tc_track([tc_track(:).orig_event_flag]), check_plot);
end

%% probabilistic tracks
gen_tracks = find(~[tc_track(:).orig_event_flag]);
if isempty(gen_tracks)
    fprintf('Input tc tracks are all historical. Please rerun with probabilistic tracks, too.\n')
    return
end

%% plot with original wind speeds (without decay)
% % control figure
if check_plot
    climada_figuresize(0.5,0.8);
    %plot([0 0],[0 150],':k')
    hold on
    xlim([-5 1000])
    ylim([900 1020])
    ylabel('Central Pressure (hPa)')
    xlabel('Distance over land (km)')
    title('Probabilistic tracks only - before pressure decay')
  
    for t_i = gen_tracks %  7:no_gen:length(tc_track)
        land_index_ = find(diff(tc_track(t_i).onLand) == 1)+1;
        sea_index_  = find(diff(tc_track(t_i).onLand) ==-1)+1;
        sea_index_  = [sea_index_ size(tc_track(t_i).onLand,2)];

        if ~isempty(land_index_)
            if length(sea_index_)<= length(land_index_)   
                % time over land
                onland_time = sea_index_ - land_index_(1:length(sea_index_));
%                onland_time(end)=onland_time(end)+1;
                for lf_i = 1:length(onland_time)
                    v_landfall  = tc_track(t_i).MaxSustainedWind(land_index_(lf_i)-1);
                    onland_dist = tc_track(t_i).distOnLand_km(land_index_(lf_i):sea_index_(lf_i)-1);
                    scale_index = find(v_landfall < v_scale_kn);
                    if ~isempty(scale_index)
                        scale_index = scale_index(1);
                        a           = onland_time(lf_i);                   
                        plot(onland_dist, tc_track(t_i).CentralPressure(land_index_(lf_i)+[0:a-1]),'.','color',cmap(scale_index,:))
                        %S_cell{scale_index}(1:a+1,end+1) = tc_track(t_i).CentralPressure(end)/p_landfall;
                    end
                end %lf_i
            end
        end
    end % loop over probabilistic tracks
end % check_plot


count_penv_lt_pcen = 0;    
for t_i = gen_tracks
    % copy before changing wind speeds
    tc_track(t_i).CentralPressure_ori = tc_track_ori_1h(t_i).CentralPressure;
    land_index_ = find(diff(tc_track(t_i).onLand) == 1)+1;
    sea_index_  = find(diff(tc_track(t_i).onLand) ==-1)+1;
    sea_index_  = [sea_index_ size(tc_track(t_i).onLand,2)];
    if ~isempty(land_index_)
        if length(sea_index_)<= length(land_index_)   
            % time over land
            onland_time = sea_index_ - land_index_(1:length(sea_index_));
            for lf_i = 1:length(onland_time)
                v_landfall  = tc_track(t_i).MaxSustainedWind(land_index_(lf_i)-1);
                p_landfall  = tc_track(t_i).CentralPressure(land_index_(lf_i)-1);
                onland_dist = tc_track(t_i).distOnLand_km(land_index_(lf_i):sea_index_(lf_i));
                scale_index = find(v_landfall < v_scale_kn);
                if ~isempty(scale_index)
                    scale_index = scale_index(1);
                    a           = onland_time(lf_i);
                    if a>1
                        
                        S = tc_track(t_i).EnvironmentalPressure(end)/min(tc_track(t_i).EnvironmentalPressure(end),p_landfall);
                        %S = 1010/min(1010,p_landfall);
                        decay  = S-(S-1).*exp(-p_rel(scale_index,1).*onland_dist);
                        tc_track(t_i).CentralPressure(land_index_(lf_i)+[0:a]) = ...
                                    p_landfall*decay;
                                    % tc_track(t_i).CentralPressure(land_index_(lf_i)+[0:a]).*decay; %%%
                                    
                        % for on-sea passages between too stretches over
                        % land, pressure is increased a bit, too account
                        % for weakening of first landfall (only affects TC
                        % on sea, not on land):
                        p_random = 0.1*(abs(randn(1)*5)+6); %mean value 1                       
                        p_diff   = diff(tc_track(t_i).CentralPressure(sea_index_(lf_i)+[-1:0]))...
                                   + p_random ;                                 
                        if sea_index_(lf_i) < length(tc_track(t_i).CentralPressure)  
                            tc_track(t_i).CentralPressure(sea_index_(lf_i):land_index_(lf_i+1)) = ...
                                            tc_track(t_i).CentralPressure(sea_index_(lf_i):land_index_(lf_i+1))-p_diff;        
%                         else
%                             tc_track(t_i).CentralPressure(sea_index_(lf_i)) = ...
%                                             tc_track(t_i).CentralPressure(sea_index_(lf_i))-p_diff; 
                        end
                    end
                end
            end %lf_i
            % clean up too small values
            % control figure
            check_plot_ = 0;% check_plot;
            if check_plot_
                fig = climada_figuresize(0.5,0.8);
                ylabel('Pressure (hPa)')
                xlabel('Time (h)')
                hold on
                for lf_i = 1:length(onland_time)
                    %plot(land_index_(lf_i)*ones(2,1),[0 150],'or-')
                    g(1) = fill([land_index_(lf_i) sea_index_(lf_i) sea_index_(lf_i) land_index_(lf_i)]-1,...
                                [0 0 150 150],[255 250 205 ]/255,'edgecolor','none');
                end
                h(1) = plot(tc_track(t_i).CentralPressure,'-');
                hold on
                h(2) = plot(tc_track(t_i).CentralPressure_ori,':k');
                legend([h g],'probabilistic track with pressure decay','probabilistic track without pressure decay','landfall',...
                             'location','southwest')
                pause
                close(fig)
            end %check_plot
        end %~isempty(land_index_)
    end
    if min(tc_track(t_i).EnvironmentalPressure-tc_track(t_i).CentralPressure)<0
        i_negative = find(tc_track(t_i).EnvironmentalPressure-tc_track(t_i).CentralPressure<0);
        tc_track(t_i).CentralPressure(i_negative)=tc_track(t_i).EnvironmentalPressure(i_negative);
        count_penv_lt_pcen = count_penv_lt_pcen+1;
    end
end %t_i
if count_penv_lt_pcen>0
    disp(['Central pressure was corrected for being larger than environmental pressure for ' num2str(count_penv_lt_pcen) ' tracks.']);
end
clear count_penv_lt_pcen i_neg* t_i onland_dist p_landfall v_landfall S decay 

if check_plot
    climada_figuresize(0.5,0.8);
    %plot([0 0],[0 150],':k')
    hold on
    xlim([-5 1000])
    ylim([900 1020])
    ylabel('Central Pressure (hPa)')
    xlabel('Distance over land (km)')
    title('Probabilistic tracks only - before pressure decay')
    for t_i = gen_tracks %  7:no_gen:length(tc_track)
        land_index_ = find(diff(tc_track(t_i).onLand) == 1)+1;
        sea_index_  = find(diff(tc_track(t_i).onLand) ==-1)+1;
        sea_index_  = [sea_index_ size(tc_track(t_i).onLand,2)];

        if ~isempty(land_index_)
            if length(sea_index_)<= length(land_index_)   
                % time over land
                onland_time = sea_index_ - land_index_(1:length(sea_index_));
%                onland_time(end)=onland_time(end)+1;
                for lf_i = 1:length(onland_time)
                    v_landfall  = tc_track(t_i).MaxSustainedWind(land_index_(lf_i)-1);
                    onland_dist = tc_track(t_i).distOnLand_km(land_index_(lf_i):sea_index_(lf_i)-1);
                    scale_index = find(v_landfall < v_scale_kn);
                    if ~isempty(scale_index)
                        scale_index = scale_index(1);
                        a           = onland_time(lf_i);                   
                        plot(onland_dist, tc_track(t_i).CentralPressure(land_index_(lf_i)+[0:a-1]),'.','color',cmap(scale_index,:))
                        %S_cell{scale_index}(1:a+1,end+1) = tc_track(t_i).CentralPressure(end)/p_landfall;
                    end
                end %lf_i
            end
        end
    end % loop over probabilistic tracks
end % check_plot


% just for control
check_plot_ = 0;%check_plot;
if check_plot_
    fig = climada_figuresize(0.5,0.8);
    for t_i = 1:no_gen:length(tc_track)
        ylabel('Central Pressure (hPa)')
        xlabel('Time (h)')
        hold on
        titlestr = sprintf('Historical track %d and %d probabilistic tracks',t_i,no_gen-1);
        title(titlestr)
        for t_ii = t_i+1:t_i+no_gen-1
            fprintf('t_ii %d \t',t_ii')
            onLand = logical(tc_track(t_ii).onLand);
            if any(~onLand)
                plot(find(~onLand), tc_track(t_ii).CentralPressure(~onLand),'b.','markersize',3)
            end
            if any(onLand)
                plot(find(onLand), tc_track(t_ii).CentralPressure(onLand),'b.')
            end
            lf = find(diff(onLand) == 1);
            fprintf('lf: ')
            fprintf(' %d, ',lf')
            fprintf('\n')
            if ~isempty(lf)
                plot(lf, tc_track(t_i).CentralPressure(lf),'bo','markersize',5)
            end
        end %t_ii
        onLand = logical(tc_track(t_i).onLand);
        if any( ~onLand)
            g(1) = plot(find(~onLand), tc_track(t_i).CentralPressure(~onLand),'k.','markersize',3);
        end
        if any(onLand)
            g(2) = plot(find(onLand), tc_track(t_i).CentralPressure(onLand),'k.');
        end
        lf   = find(diff(onLand) == 1);
        if ~isempty(lf)
            g(3) = plot(lf, tc_track(t_i).CentralPressure(lf),'or','markersize',7);
        end
        ylim([0 max(tc_track(t_i).CentralPressure)*1.2])
        xlim([0 size(onLand,2)*1.1])
        legend(g,'wind on sea', 'wind on land', 'landfall','location','southwest')
        pause
        delete(gca)
        clear g
    end %t_i
end







