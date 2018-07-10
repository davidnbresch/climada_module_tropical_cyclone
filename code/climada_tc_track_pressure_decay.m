function tc_track = climada_tc_track_pressure_decay(tc_track, p_rel, check_plot)

% NAME:
%   climada_tc_track_pressure_decay
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
%          %   pressure decay = S-(S-1)*exp(-A*x), where A = p_rel(1,1), S=LastPressureOfTrack/PressureAtLandfall=p_rel(1,2)
%          A = p_rel(1,1), S = p_rel(1,2), can be newly calculated or can
%          be loaded from data within globalGDP modul
%   check_plot: to create plot
% OUTPUTS:
%   same structure now CentralPressure in tc track decays after landfall
% RESTRICTIONS:
%   changes CentralPressure only.
% MODIFICATION HISTORY:
% Samuel Eberenz, eberenz@posteo.eu, 20180709, adapted from climada_tc_track_wibd_decay_calculate to change CentralPressure p instead of MaxSustainedWind v (not yet fully documented)
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

% refine tc_tracks to 1 h
tc_track = climada_tc_equal_timestep(tc_track,1);


% make a copy
tc_track_ori_1h = tc_track;



%% find nodes on land and over sea
tc_track = climada_tc_on_land(tc_track);


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
cmap         = jet(no_cat);

%% calculate exponential decay of wind speed after landfall if not given
if isempty(p_rel)
    fprintf('Calculate p rel (parameters for exponential wind decay) based on historical tracks...\n')
    [~, p_rel] = climada_tc_track_pressure_decay_calculate(tc_track([tc_track(:).orig_event_flag]), check_plot);
end

if check_plot
    % create figure with fitted functions, relative decay
    h = zeros(1,no_cat);
    climada_figuresize(0.5,0.8);
    hold on
    xlim([-5 150])
    ylim([0.9 1.3])
    ylabel('Relative pressure (on landfall = 1)')
    timestep  = datenum(0,0,diff(tc_track(1).datenum(1:2)))*24;
    xlabelstr = sprintf('Time steps after landfall (h)');
    xlabel(xlabelstr)
    cmap   = jet(no_cat);
    for cat_i = 1:no_cat
        if ~isempty(p_vs_lf_time{cat_i})
            hg = plot( [0:size(p_vs_lf_time{cat_i},1)-1],...
                              v_vs_lf_time_relative{cat_i},'.','color',cmap(cat_i,:),'markersize',5);
            h(cat_i) = hg(1);              
        end
    end
    %p_rel = p_rel*2;
    x_fit = [1:0.5:150];
    for cat_i = 1:no_cat
        y_fit = polyval(p_rel(cat_i,:),x_fit);
        g = plot(x_fit,exp(y_fit),'-','color',cmap(cat_i,:));
        h(cat_i) = g(1);
        %semilogy(v_vs_lf_time{cat_i},'.','color',cmap(cat_i,:))
        hold on
    end    
    legendstr = {'Tropical Depression','Tropical Storm','Hurrican Cat. 1','Hurrican Cat. 2','Hurrican Cat. 3','Hurrican Cat. 4','Hurrican Cat. 5'};
    for cat_i = 1:no_cat
        legendstr{cat_i} = sprintf('%s,   y = %1.3f-%1.3f*exp(%1.4f x)',legendstr{cat_i},p_rel(cat_i,2),p_rel(cat_i,2)-1, -p_rel(cat_i,1));
    end
    legend(h,legendstr)
    title('Relative pressure decay in relation to time after landfall')
end %check_plot


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
    plot([0 0],[0 150],':k')
    hold on
    xlim([-5 150])
    ylim([900 1020])
    ylabel('Central Pressure (hPa)')
    xlabel('Time after landfall (h)')
    title('Probabilistic tracks only - before pressure decay')
    for t_i = gen_tracks %  7:no_gen:length(tc_track)
        land_index_ = find(diff(tc_track(t_i).onLand) == 1)+1;
        sea_index_  = find(diff(tc_track(t_i).onLand) ==-1)+1;
        sea_index_  = [sea_index_ size(tc_track(t_i).onLand,2)];

        if ~isempty(land_index_)
            if length(sea_index_)<= length(land_index_)   
                % time over land
                onland_time = sea_index_ - land_index_(1:length(sea_index_));

                for lf_i = 1:length(onland_time)
                    v_landfall  = tc_track(t_i).MaxSustainedWind(land_index_(lf_i)-1);
                    
                    scale_index = find(v_landfall < v_scale_kn);
                    if ~isempty(scale_index)
                        scale_index = scale_index(1);
                        a           = onland_time(lf_i);                   
                        plot(0:a, tc_track(t_i).CentralPressure(land_index_(lf_i)+[0:a]),'.','color',cmap(scale_index,:))
                        %S_cell{scale_index}(1:a+1,end+1) = tc_track(t_i).CentralPressure(end)/p_landfall;
                    end
                end %lf_i
            end
        end
    end % loop over probabilistic tracks
end % check_plot


    
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
                scale_index = find(v_landfall < v_scale_kn);
                if ~isempty(scale_index)
                    scale_index = scale_index(1);
                    a           = onland_time(lf_i);
                    if a>1
                        %decay    = exp(polyval(p_rel(scale_index,:),1:a));
                        decay  = p_rel(scale_index,2)-(p_rel(scale_index,2)-1).*exp(-p_rel(scale_index,1).*[1:a]);
                        tc_track(t_i).CentralPressure(land_index_(lf_i)-1+[1:a]) = ...
                                                             p_landfall*decay;
                        p_random = 0.1*(abs(randn(1)*5)+6); %mean value 1                       
                        p_diff   = diff(tc_track(t_i).CentralPressure(sea_index_(lf_i)+[-1:0]))...
                                   - p_random ;                                 
                        if sea_index_(lf_i) < length(tc_track(t_i).CentralPressure)  
                            tc_track(t_i).CentralPressure(sea_index_(lf_i):land_index_(lf_i+1)) = ...
                                            tc_track(t_i).CentralPressure(sea_index_(lf_i):land_index_(lf_i+1))-p_diff;        
                        else
                            tc_track(t_i).CentralPressure(sea_index_(lf_i)) = ...
                                            tc_track(t_i).CentralPressure(sea_index_(lf_i))-p_diff; 
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
end %t_i

% no numbers above 1015
for t_i = 1:length(tc_track)
    if any(tc_track(t_i).CentralPressure>1015)
        tc_track(t_i).CentralPressure(tc_track(t_i).CentralPressure>1015) = 1015;
    end
end


if check_plot
    % % control figure
    climada_figuresize(0.5,0.8);
    plot([0 0],[0 150],':k')
    hold on
    xlim([-5 150])
    ylim([900 1020])
    ylabel('Pressure (hPa)')
    xlabel('Time after landfall (h)')
    title('Probabilistic tracks only - after pressure decay at landfall')
    for t_i = gen_tracks %  7:no_gen:length(tc_track)
        land_index_ = find(diff(tc_track(t_i).onLand) == 1)+1;
        sea_index_  = find(diff(tc_track(t_i).onLand) ==-1)+1;
        sea_index_  = [sea_index_ size(tc_track(t_i).onLand,2)];

        if ~isempty(land_index_)
            if length(sea_index_)<= length(land_index_)   
                % time over land
                onland_time = sea_index_ - land_index_(1:length(sea_index_));

                for lf_i = 1:length(onland_time)
                    v_landfall  = tc_track(t_i).MaxSustainedWind(land_index_(lf_i)-1);
                    %p_landfall  = tc_track(t_i).CentralPressure(land_index_(lf_i)-1);
                    scale_index = find(v_landfall < v_scale_kn);
                    if ~isempty(scale_index)
                        scale_index = scale_index(1);
                        a           = onland_time(lf_i);                   
                        plot(0:a, tc_track(t_i).CentralPressure(land_index_(lf_i)+[0:a]),'.','color',cmap(scale_index,:))
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







