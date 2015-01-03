function tc_track = climada_tc_wind_decay(tc_track)
% UNDOCUMENTED
%-

% init global variables
global climada_global
if ~climada_init_vars, return; end

% check inputs, and set default values
if ~exist('tc_track'       , 'var'), tc_track      = []  ; end

% prompt for tc_track if not given
if isempty(tc_track)
    tc_track             = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    tc_track_default     = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select HISTORICAL tc track .mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select HISTORICAL tc track:',tc_track_default);
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

% make a copy
tc_track_ori = tc_track;

%% refine tc_tracks to 1 h
msgstr   = sprintf('refining %i tracks to 1h timestep\n',length(tc_track));
if climada_global.waitbar,h        = waitbar(0,msgstr);end
mod_step = 100;
for t_i = 1:length(tc_track_ori)
    tc_track(t_i) = climada_tc_equal_timestep(tc_track_ori(t_i),1);
    if mod(t_i,mod_step)==0 && climada_global.waitbar
        mod_step          = 100;
        msgstr            = sprintf('Refining tracks to 1h timestep\n%i/%i tracks',t_i, length(tc_track));
        waitbar(t_i/length(tc_track),h,msgstr); % update waitbar
    end
end
if climada_global.waitbar,close(h);end
tc_track_ori_1h = tc_track;



%% find nodes on land and over sea
% set modul data directory
modul_data_dir = [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
load([modul_data_dir filesep 'border_mask_10km'])

msgstr   = sprintf('processing %i tracks\n',length(tc_track));
if climada_global.waitbar,h        = waitbar(0,msgstr);end
mod_step = 100;
for t_i = 1:length(tc_track)
    i = round((tc_track(t_i).lon - (border_mask.lon_range(1)-border_mask.resolution_x/2)) /border_mask.resolution_x)+1;
    j = round((tc_track(t_i).lat - (border_mask.lat_range(1)-border_mask.resolution_y/2)) /border_mask.resolution_y)+1;
    i(i>size(border_mask.world_mask,2)) = size(border_mask.world_mask,2);
    j(j>size(border_mask.world_mask,1)) = size(border_mask.world_mask,1);
    for n_i = 1:length(i)
        tc_track(t_i).onLand(n_i) = border_mask.world_mask(j(n_i),i(n_i));
    end 
    if mod(t_i,mod_step)==0 && climada_global.waitbar
        mod_step          = 100;
        msgstr            = sprintf('Add onLand variable to each tc_track \n%i/%i tracks',t_i, length(tc_track));
        waitbar(t_i/length(tc_track),h,msgstr); % update waitbar
    end
end
if climada_global.waitbar,close(h);end

% number of generated and historical tracks
no_ori = sum([tc_track(:).orig_event_flag]);
no_gen = length(tc_track)/no_ori;


%% plot wind speed after land fall of historical tracks
% climada_figuresize(0.5,0.8);
% plot([0 0],[0 150],':k')
% hold on
% xlim([0 150])
% ylim([0 150])
% ylabel('Wind speed (kn)')
% xlabel('Time after landfall (h)')
v_scale_kn   = [34 64 83 96 113 135 500];
no_cat       = size(v_scale_kn,2);
v_vs_lf_time = cell(1,no_cat);

for t_i = 1:no_gen:length(tc_track)
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
                    v_vs_lf_time{scale_index}(1:a+1,end+1) = tc_track(t_i).MaxSustainedWind(land_index_(lf_i)-1+[0:a])';
                    %plot(1:a, tc_track(t_i).MaxSustainedWind(land_index+[0:a-1]),'-')
                end
            end %lf_i
        end
    end
end   
% put into one structure and set zeros to nan
for cat_i = 1:no_cat
    v_vs_lf_time{cat_i}(v_vs_lf_time{cat_i} == 0) = nan;
end


% % create figure
% climada_figuresize(0.5,0.8);
% % plot([0 0],[0 150],':k')
% hold on
% xlim([-5 150])
% ylim([0 150])
% ylabel('Wind speed (kn)')
% xlabel('Time after landfall (h)')
% no_cat = size(v_vs_lf_time,2);
% cmap   = jet(no_cat);
% for cat_i = 1:no_cat
%     g = plot( [0:size(v_vs_lf_time{cat_i},1)-1],...
%               v_vs_lf_time{cat_i},'.','color',cmap(cat_i,:));
%     h(cat_i) = g(1);
%     %semilogy(v_vs_lf_time{cat_i},'.','color',cmap(cat_i,:))
%     hold on
% end
% legend(h,'Tropical Depression','Tropical Storm','Hurrican Cat. 1','Hurrican Cat. 2','Hurrican Cat. 3','Hurrican Cat. 4','Hurrican Cat. 5' )

% % fit exponential decay
% p = zeros(no_cat,2);
% for cat_i = 1:no_cat
%     x = repmat([1:size(v_vs_lf_time{cat_i},1)]',size(v_vs_lf_time{cat_i},2),1);
%     y = v_vs_lf_time{cat_i}(:);
%     p(cat_i,:) = polyfit(x(~isnan(y)), log(y(~isnan(y))),1);
% end

% % create figure with fitted functions
% climada_figuresize(0.5,0.8);
% hold on
% xlim([-5 150])
% ylim([0 150])
% ylabel('Wind speed (kn)')
% xlabel('Time after landfall (h)')
% no_cat = size(v_vs_lf_time,2);
% cmap   = jet(no_cat);
% for cat_i = 1:no_cat
%     plot( [0:size(v_vs_lf_time{cat_i},1)-1],...
%            v_vs_lf_time{cat_i},'.','color',cmap(cat_i,:),'markersize',5);
% end
% for cat_i = 1:no_cat
%     x_fit = [0:2:150];
%     y_fit = polyval(p(cat_i,:),x_fit);
%     g = plot(x_fit,exp(y_fit),'.-','color',cmap(cat_i,:));
%     h(cat_i) = g(1);
%     %semilogy(v_vs_lf_time{cat_i},'.','color',cmap(cat_i,:))
%     hold on
% end
% legend(h,'Tropical Depression','Tropical Storm','Hurrican Cat. 1','Hurrican Cat. 2','Hurrican Cat. 3','Hurrican Cat. 4','Hurrican Cat. 5' )
% title('Wind speed decay in relation to time after landfall')


%% relative decay for different categories
for cat_i = 1:no_cat
    v_vs_lf_time_relative{cat_i} = bsxfun(@rdivide,v_vs_lf_time{cat_i},v_vs_lf_time{cat_i}(1,:));
end

%% fit exponential decay with intercept 1 at landfall
p_rel = zeros(no_cat,2);
for cat_i = 1:no_cat
    x = repmat([1:size(v_vs_lf_time_relative{cat_i},1)]',size(v_vs_lf_time_relative{cat_i},2),1);
    y = v_vs_lf_time_relative{cat_i}(:);
    B = log(y(~isnan(y))) ./ x(~isnan(y)) ;
    p_rel(cat_i,1) = mean(B);
end
nan_index = isnan(p_rel(:,1));
n_nan     = find(~nan_index);
for cat_i = 1:no_cat
    if isnan(p_rel(cat_i,1))
        [c closest]    = min(abs(cat_i-n_nan));
        p_rel(cat_i,1) = p_rel(n_nan(closest(1)),1);
        fprintf('No historical track Category %d. Take wind decay parameters from Category %d\n', cat_i-2, n_nan(closest(1))-2)
    end
end


% create figure with fitted functions, relative decay
climada_figuresize(0.5,0.8);
hold on
xlim([-5 150])
ylim([0 1.1])
ylabel('Relative wind speed (on landfall = 1)')
timestep  = datenum(0,0,diff(tc_track(1).datenum(1:2)))*24;
xlabelstr = sprintf('Time after landfall (h)');
xlabel(xlabelstr)
cmap   = jet(no_cat);
for cat_i = 1:no_cat
    if ~isempty(v_vs_lf_time{cat_i})
        plot( [0:size(v_vs_lf_time{cat_i},1)-1],...
              v_vs_lf_time_relative{cat_i},'.','color',cmap(cat_i,:),'markersize',5);
    end
end
for cat_i = 1:no_cat
    x_fit = [1:0.5:150];
    y_fit = polyval(p_rel(cat_i,:),x_fit);
    g = plot(x_fit,exp(y_fit),'-','color',cmap(cat_i,:));
    h(cat_i) = g(1);
    %semilogy(v_vs_lf_time{cat_i},'.','color',cmap(cat_i,:))
    hold on
end
legend(h,'Tropical Depression','Tropical Storm','Hurrican Cat. 1','Hurrican Cat. 2','Hurrican Cat. 3','Hurrican Cat. 4','Hurrican Cat. 5' )
title('Relative wind speed decay in relation to time after landfall')



%% probabilistic tracks - wind speed decay
v_threshold = 0;
gen_tracks = find(~[tc_track(:).orig_event_flag]);
for t_i = gen_tracks
    % copy before changing wind speeds
    tc_track(t_i).MaxSustainedWind_ori = tc_track_ori(t_i).MaxSustainedWind;
    
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
                        decay       = exp(polyval(p_rel(cat_i,:),1:a));
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
%             fig = climada_figuresize(0.5,0.8);
%             hold on
%             for lf_i = 1:length(onland_time)
%                 %plot(land_index_(lf_i)*ones(2,1),[0 150],'or-')
%                 fill([land_index_(lf_i) sea_index_(lf_i) sea_index_(lf_i) land_index_(lf_i)]-1,...
%                      [0 0 150 150],[255 250 205 ]/255,'edgecolor','none')
%             end
%             plot(tc_track(t_i).MaxSustainedWind,'-')
%             hold on
%             plot(tc_track(t_i).MaxSustainedWind_ori,':k')
%             pause
%             close(fig)
        end
    end
end %t_i



% % control figure
climada_figuresize(0.5,0.8);
plot([0 0],[0 150],':k')
hold on
xlim([-5 150])
ylim([0 150])
ylabel('Wind speed (kn)')
xlabel('Time after landfall (h)')
title('Probabilistic tracks only')

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
end   




