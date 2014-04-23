function tc_track = climada_tc_track_wind_decay(tc_track, p_rel, check_plot)

% Incorporate wind decay after landfall for probabilistic tracks, based on
% historical tc tracks
% NAME:
%   climada_tc_track_wind_decay
% PURPOSE:
%   incorporate wind decay after landfall of probabilistic tracks. The
%   decay can be calculated based on historical tracks or can be loaded
%   from a mat file (p_rel)
%   within:  climada_tc_random_walk_position_windspeed
% CALLING SEQUENCE:
%   tc_track = climada_tc_track_wind_decay(tc_track, p_rel, check_plot)
% EXAMPLE:
%   tc_track = climada_tc_track_wind_decay(tc_track)
% INPUTS:
%   none, if tc_track empty prompted for
% OPTIONAL INPUT PARAMETERS:
%   p_rel: parameters for exponential decay, y = exp(B) * exp(x*A), where 
%          A = p_rel(1,1), B = p_rel(1,2), can be newly calculated or can
%          be loaded from data within globalGDP modul
%   check_plot: to create plot
% OUTPUTS:
%   same structure now MaxSustainedWind in tc track decays after landfall
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, 20121203
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


%% refine tc_tracks to 1 h
msgstr   = sprintf('refining %i tracks to 1h timestep\n',length(tc_track));
h        = waitbar(0,msgstr);
mod_step = 500;
for t_i = 1:length(tc_track)
    tc_track(t_i) = climada_tc_equal_timestep(tc_track(t_i),1);
    if mod(t_i,mod_step)==0
        mod_step          = 500;
        msgstr            = sprintf('Refining tracks to 1h timestep\n%i/%i tracks',t_i, length(tc_track));
        waitbar(t_i/length(tc_track),h,msgstr); % update waitbar
    end
end
close(h)
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
v_vs_lf_time = cell(1,no_cat);
cmap         = jet(no_cat);

%% calculate exponential decay of wind speed after landfall if not given
if isempty(p_rel)
    fprintf('Calculate p rel (parameters for exponential wind decay) based on historical tracks...\n')
    [unused p_rel] = climada_tc_track_wind_decay_calculate(tc_track([tc_track(:).orig_event_flag]), check_plot);
end

if check_plot
    % create figure with fitted functions, relative decay
    h = zeros(1,no_cat);
    climada_figuresize(0.5,0.8);
    hold on
    xlim([-5 150])
    ylim([0 1.1])
    ylabel('Relative wind speed (on landfall = 1)')
    timestep  = datenum(0,0,diff(tc_track(1).nodetime_mat(1:2)))*24;
    xlabelstr = sprintf('Time steps after landfall (h)');
    xlabel(xlabelstr)
    cmap   = jet(no_cat);
    for cat_i = 1:no_cat
        if ~isempty(v_vs_lf_time{cat_i})
            hg = plot( [0:size(v_vs_lf_time{cat_i},1)-1],...
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
        legendstr{cat_i} = sprintf('%s, y = exp(%1.4f t)',legendstr{cat_i}, p_rel(cat_i,1));
    end
    legend(h,legendstr)
    title('Relative wind speed decay in relation to time after landfall')
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
    ylim([0 150])
    ylabel('Wind speed (kn)')
    xlabel('Time after landfall (h)')
    title('Probabilistic tracks only - before wind speed decay')
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
                        plot(0:a, tc_track(t_i).MaxSustainedWind(land_index_(lf_i)+[0:a]),'.','color',cmap(scale_index,:))
                    end
                end %lf_i
            end
        end
    end % loop over probabilistic tracks
end % check_plot


    
for t_i = gen_tracks
    % copy before changing wind speeds
    tc_track(t_i).MaxSustainedWind_ori = tc_track_ori_1h(t_i).MaxSustainedWind;
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
                    if a>1
                        decay    = exp(polyval(p_rel(scale_index,:),1:a));
                        tc_track(t_i).MaxSustainedWind(land_index_(lf_i)-1+[1:a]) = ...
                                                             v_landfall*decay;
                        v_random = abs(randn(1)*5)+6; %mean value 10                       
                        v_diff   = diff(tc_track(t_i).MaxSustainedWind(sea_index_(lf_i)+[-1:0]))...
                                   - v_random ;                                 
                        if sea_index_(lf_i) < length(tc_track(t_i).MaxSustainedWind)  
                            tc_track(t_i).MaxSustainedWind(sea_index_(lf_i):land_index_(lf_i+1)) = ...
                                            tc_track(t_i).MaxSustainedWind(sea_index_(lf_i):land_index_(lf_i+1))-v_diff;        
                        else
                            tc_track(t_i).MaxSustainedWind(sea_index_(lf_i)) = ...
                                            tc_track(t_i).MaxSustainedWind(sea_index_(lf_i))-v_diff; 
                        end
                    end
                end
            end %lf_i
            % clean up too small values
            % control figure
            check_plot_ = 0;
            if check_plot_
                fig = climada_figuresize(0.5,0.8);
                ylabel('Wind speed (kn)')
                xlabel('Time (h)')
                hold on
                for lf_i = 1:length(onland_time)
                    %plot(land_index_(lf_i)*ones(2,1),[0 150],'or-')
                    g(1) = fill([land_index_(lf_i) sea_index_(lf_i) sea_index_(lf_i) land_index_(lf_i)]-1,...
                                [0 0 150 150],[255 250 205 ]/255,'edgecolor','none');
                end
                h(1) = plot(tc_track(t_i).MaxSustainedWind,'-');
                hold on
                h(2) = plot(tc_track(t_i).MaxSustainedWind_ori,':k');
                legend([h g],'probabilistic track with wind decay','probabilistic track without wind decay','landfall',...
                             'location','southwest')
                pause
                close(fig)
            end %check_plot
        end %~isempty(land_index_)
    end
end %t_i

% no negativ numbers
for t_i = 1:length(tc_track)
    if any(tc_track(t_i).MaxSustainedWind<0)
        tc_track(t_i).MaxSustainedWind(tc_track(t_i).MaxSustainedWind<0) = 0;
    end
end


if check_plot
    % % control figure
    climada_figuresize(0.5,0.8);
    plot([0 0],[0 150],':k')
    hold on
    xlim([-5 150])
    ylim([0 150])
    ylabel('Wind speed (kn)')
    xlabel('Time after landfall (h)')
    title('Probabilistic tracks only - after wind decay at landfall')
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
                        plot(0:a, tc_track(t_i).MaxSustainedWind(land_index_(lf_i)+[0:a]),'.','color',cmap(scale_index,:))
                    end
                end %lf_i
            end
        end
    end % loop over probabilistic tracks
end % check_plot


% just for control
check_plot_ = 0;
if check_plot_
    fig = climada_figuresize(0.5,0.8);
    for t_i = 1:no_gen:length(tc_track)
        ylabel('Wind speed (kn)')
        xlabel('Time (h)')
        hold on
        titlestr = sprintf('Historical track %d and %d probabilistic tracks',t_i,no_gen-1);
        title(titlestr)
        for t_ii = t_i+1:t_i+no_gen-1
            fprintf('t_ii %d \t',t_ii')
            onLand = logical(tc_track(t_ii).onLand);
            if any(~onLand)
                plot(find(~onLand), tc_track(t_ii).MaxSustainedWind(~onLand),'b.','markersize',3)
            end
            if any(onLand)
                plot(find(onLand), tc_track(t_ii).MaxSustainedWind(onLand),'b.')
            end
            lf = find(diff(onLand) == 1);
            fprintf('lf: ')
            fprintf(' %d, ',lf')
            fprintf('\n')
            if ~isempty(lf)
                plot(lf, tc_track(t_i).MaxSustainedWind(lf),'bo','markersize',5)
            end
        end %t_ii
        onLand = logical(tc_track(t_i).onLand);
        if any( ~onLand)
            g(1) = plot(find(~onLand), tc_track(t_i).MaxSustainedWind(~onLand),'k.','markersize',3);
        end
        if any(onLand)
            g(2) = plot(find(onLand), tc_track(t_i).MaxSustainedWind(onLand),'k.');
        end
        lf   = find(diff(onLand) == 1);
        if ~isempty(lf)
            g(3) = plot(lf, tc_track(t_i).MaxSustainedWind(lf),'or','markersize',7);
        end
        ylim([0 max(tc_track(t_i).MaxSustainedWind)*1.2])
        xlim([0 size(onLand,2)*1.1])
        legend(g,'wind on sea', 'wind on land', 'landfall','location','southwest')
        pause
        delete(gca)
        clear g
    end %t_i
end







