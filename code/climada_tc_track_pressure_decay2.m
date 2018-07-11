function tc_track = climada_tc_track_pressure_decay2(tc_track, timelim_h, check_plot)

% NAME:
%   climada_tc_track_pressure_decay2
% MODULE:
%   tropical_cyclone
% PURPOSE:
%   Incorporate strong pressure decay after landfall for probabilistic tracks, based 
%   purely on time after landfall. 
%   After landfall, central pressure is increased exponentially 
%   to reach 1015mb N hours after land fall (ususally N=48)
%   for:  isimip_windfield_holland
% CALLING SEQUENCE:
%   tc_track = climada_tc_track_pressure_decay2(tc_track, timelim_h, check_plot)
%   hazard   = isimip_tc_hazard_set(tc_track_decay,hazard_file,centroids);
% EXAMPLE:
%   tc_track = climada_tc_track_pressure_decay2([], 72, 1)
% INPUTS:
%   none, if tc_track empty prompted for
% OPTIONAL INPUT PARAMETERS:
%   timelim_h: time in hours after landfall when pressure should reach environmental pressure.
%           (default = 48)
%          pressure after landfall is computed from each track and
%          timelim_h to get an exponential saturation:
%   check_plot: to create plot
% OUTPUTS:
%   same structure now CentralPressure in tc track decays after landfall
% RESTRICTIONS:
%   changes CentralPressure only.
% MODIFICATION HISTORY:
% Samuel Eberenz, eberenz@posteo.eu, 20180711, adapted from climada_tc_track_wind_decay to change CentralPressure p 
%-

% init global variables
global climada_global
if ~climada_init_vars, return; end

% check inputs, and set default values
if ~exist('tc_track'       , 'var'), tc_track      = []  ; end
%if ~exist('p_rel'          , 'var'), p_rel         = []  ; end
if ~exist('timelim_h'       , 'var'), timelim_h      = []  ; end
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
if isempty(timelim_h),timelim_h=48;end % hours after landfall
% refine tc_tracks to 1 h
tc_track = climada_tc_equal_timestep(tc_track,1);

% make a copy
tc_track_ori_1h = tc_track;


%% find nodes on land and over sea
tc_track = climada_tc_on_land(tc_track);

% number of generated and historical tracks
no_ori = sum([tc_track(:).orig_event_flag]);
no_gen = length(tc_track)/no_ori;

fsp = 1020; % mb, fixed_saturation_pressure
mtp = 1015; % mb, max_target_pressure 


v_scale_kn   = [34 64 83 96 113 135 500];
no_cat       = size(v_scale_kn,2);
p_vs_lf_time = cell(1,no_cat);
cmap         = jet(no_cat);


%% probabilistic tracks
% gen_tracks = find(~[tc_track(:).orig_event_flag]);
% if isempty(gen_tracks)
%     fprintf('Input tc tracks are all historical. Please rerun with probabilistic tracks, too.\n')
%     return
% end

%% plot with original wind speeds (without decay)
% % control figure
if check_plot
    climada_figuresize(0.5,0.8);
    plot([0 0],[0 min(94,timelim_h+12)],':k')
    hold on
    xlim([-5 min(94,timelim_h+12)])
    ylim([900 1020])
    ylabel('Central Pressure (hPa)')
    xlabel('Time after landfall (h)')
    title('All tracks - before pressure decay')
    for t_i = 1:length(tc_track) %gen_tracks % 
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


    
for t_i = 1:length(tc_track)
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
                        % fsp = 1020;% mb (fixed_saturation_pressure)
                        % mtp = 1015; % (max_target_pressure)
                        % A = (-1/timelim_h) * log((fsp-min(mtp,tc_track(t_i).EnvironmentalPressure(end)))/(fsp-min(p_landfall,mtp-5)));
                        A = (-1/timelim_h) * log((fsp-mtp)/(fsp-min(p_landfall,mtp-5)));
                        % i.e.  A = -1/48 * ln((1020-p_env)/(1020-p0))
                        % decay = fsp-(fsp-p_landfall)*exp(-real(A).*[1:a]);
                        % p(t) = S-(S-p0)*exp(-A*t)
                        
                        tc_track(t_i).CentralPressure(land_index_(lf_i)-1+[1:a]) = ...
                            fsp-(fsp-p_landfall)*exp(-real(A).*[1:a]); % pressure saturating towards target pressure
                        
                    %       p_landfall*decay;
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
                    %plot(land_index_(lf_i)*ones(2,1),[0 min(94,timelim_h+12)],'or-')
                    g(1) = fill([land_index_(lf_i) sea_index_(lf_i) sea_index_(lf_i) land_index_(lf_i)]-1,...
                                [0 0 150 150],[255 250 205 ]/255,'edgecolor','none');
                end
                h(1) = plot(tc_track(t_i).CentralPressure,'-');
                hold on
                h(2) = plot(tc_track(t_i).CentralPressure_ori,':k');
                legend([h g],'track with pressure decay','track without pressure decay','landfall',...
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
    plot([0 0],[0 min(94,timelim_h+12)],':k')
    hold on
    xlim([-5 min(94,timelim_h+12)])
    ylim([900 1020])
    ylabel('Pressure (hPa)')
    xlabel('Time after landfall (h)')
    title('All tracks - pressure decay at landfall')
    for t_i = 1:length(tc_track) %  7:no_gen:length(tc_track)
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







