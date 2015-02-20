function h = climada_plot_IFC_return(hazard, centroids, important_centroid, check_printplot, check_xlog)
% climada intensity vs return period for multiple centroids
% based on probabistic hazard set
% NAME:
%   climada_plot_IFC_return
% PURPOSE:
%   find all wind speed greater zero that hit important centroid
%   (historical and probabilistic storms)
%   calculate the exceedence frequency and return period 
%   fit a gumbel distribution for requested return period and extrapolate
%   according wind speed values
%   plot wind speed vs return period
% CALLING SEQUENCE:
%   climada_plot_IFC_return(hazard,centroids,important_centroid,check_printplot);
% EXAMPLE:
%   climada_plot_IFC_return
%   climada_plot_IFC_return(hazard, centroids, [7 23 70000])
% INPUTS:
%   hazard: probabilistic hazard set (structure)
%   centroids
%   important_centroid: index number of one or multiple centroids (scalar or array)
% OPTIONAL INPUT PARAMETERS:
%   check_printplot: if set to 1 will print (save) figure
% OUTPUTS:
%   none
% RESTRICTIONS:
%   none
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20110519
% Lea Mueller, muellele@gmail.com, 20150123, changed hazard.arr to hazard.intensity
% Gilles Stassen, gillesstassen@hotmail.com, 20150123, ylabel for different perils
% Gilles Stassen, gillesstassen@hotmail.com, 20150128, fixed an issue with using 'pos' for indexing; 
%                                                      incorporate non-uniform frequency vector compatability
% Gilles Stassen, gillesstassen@hotmail.com, 20150129, added logical check_logx argument to trigger logarithmic x axis
%-
                                                               

%% 
global climada_global
if ~climada_init_vars,return;end % init/import global variables
if ~exist('hazard'            ,'var'), hazard             = []; end
if ~exist('centroids'         ,'var'), centroids          = []; end
if ~exist('important_centroid','var'), important_centroid = []; end
if ~exist('check_printplot'   ,'var'), check_printplot    = []; end
<<<<<<< HEAD
if ~exist('check_xlog'        ,'var'), check_xlog         = []; end
=======
if ~exist('check_xlog'        ,'var'), check_xlog         = 0 ; end
>>>>>>> 5987dc3e92369aef9b1219a4affe9d03cdc4b25a

h =[];

%% prompt for hazard_set if not given
if isempty(hazard) % local GUI
    hazard               = [climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    default_hazard       = [climada_global.data_dir filesep 'hazards' filesep 'select hazard .mat'];
    [filename, pathname] = uigetfile(hazard, 'Select hazard event set for EDS calculation:',default_hazard);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard=fullfile(pathname,filename);
    end
end
% load the hazard set, if a filename has been passed
if ~isstruct(hazard)
    hazard_file=hazard;hazard=[];
    load(hazard_file);
end

%% prompt for centroids if not given
if isempty(centroids)
    centroids = [climada_global.system_dir filesep '*.mat'];
    [filename, pathname] = uigetfile(centroids, 'Select centroids:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        centroids = fullfile(pathname,filename);
    end
end
% load the centroids, if a filename has been passed
if ~isstruct(centroids)
    centroids_file = centroids;
    centroids      = [];
    vars = whos('-file', centroids_file);
    load(centroids_file);
    if ~strcmp(vars.name,'centroids')
        centroids = eval(vars.name);
        clear (vars.name)
    end
end

%% prompt for important_centroid if not given
if isempty(important_centroid)
    prompt   ='Type specific No. of centroids [e.g. 1, 7, 34]:';
    name     =' No. of centroid';
    defaultanswer = {'7'};
    answer = inputdlg(prompt,name,1,defaultanswer);
    answer = cell2mat(answer);
    important_centroid = sscanf(answer,'%d');
    %important_centroid = str2num(answer{1}); 
end

indx_invalid = important_centroid>length(hazard.lon);
if any(indx_invalid)
    important_centroid(indx_invalid)=length(hazard.lon);
end

name = cell(1,length(important_centroid));
if isfield(centroids,'names')
    for imp_i = 1:length(important_centroid)
        name{imp_i} = centroids.names{important_centroid(imp_i)};
    end
else
    for imp_i = 1:length(important_centroid)
        name{imp_i} = ['Centroid ' int2str(important_centroid(imp_i))];
    end
end
titlestr = [];
for imp_i = 1:length(important_centroid)
    titlestr = [titlestr name{imp_i} ' ']; 
end

%%
threshold                = 5;
no_generated             = hazard.event_count / hazard.orig_event_count;
requested_return_periods = [1 5 10 20 50 100 250 1000];
r_length                 = 1:1:length(requested_return_periods);


for imp_i = 1:length(important_centroid)
    %% historical data
    %1: intensity
    [int{imp_i},int_ndx] = sort(full(hazard.intensity(1:no_generated:end,important_centroid(imp_i))),'descend');
    neg                = int{imp_i} < threshold;
    pos                = int{imp_i} >= threshold;
    %frequency
    int{imp_i}(:,3)    = hazard.frequency(int_ndx).*no_generated;
    %3: exceedence frequency
    int{imp_i}(:,3)    = cumsum(int{imp_i}(:,3));
    int{imp_i}(neg,:)  = nan;
    %2: fitted intensity
    p(:,imp_i)         = polyfit(log(int{imp_i}(pos,3) ), int{imp_i}(pos,1), 1); 
    int{imp_i}(:,2)    = polyval(p(:,imp_i), log(int{imp_i}(:,3)));
    neg                = int{imp_i}(:,2) < 0;
    int{imp_i}(neg,2)  = nan;
    
    %% probabilistic data
    %4: intensity
    [int_, int_ndx_]   = sort(full(hazard.intensity(:,important_centroid(imp_i))),'descend');
    int{imp_i}(:,4)    = int_(1:length(int{imp_i}));
    neg                = int{imp_i}(:,4) < threshold;
    %frequency
    int{imp_i}(:,5)    = hazard.frequency(int_ndx_(1:length(int{imp_i})));
    %5: exceedence frequency
    int{imp_i}(:,5)    = cumsum(int{imp_i}(:,5));
    int{imp_i}(neg,4:5)= nan;
    
    %% fit a Gumbel-distribution
    %8: exceedence frequency
    int{imp_i}(r_length,8)= 1./requested_return_periods;
    %6: intensity hist.
    int{imp_i}(r_length,6)= polyval(p(:,imp_i), log(int{imp_i}(r_length,8)));
    neg                   = int{imp_i}(:,6) < 0;
    int{imp_i}(neg,6)     = nan;
    %7: intensity prob.
    pos                   = int{imp_i}(:,4) >= threshold;
    p_(:,imp_i)           = polyfit(log(int{imp_i}(pos,5) ), int{imp_i}(pos,4), 1); 
    int{imp_i}(r_length,7)= polyval(p_(:,imp_i), log(int{imp_i}(r_length,8)));
    neg                   = int{imp_i}(:,7) < 0;
    int{imp_i}(neg,7)     = nan;
    int{imp_i}(r_length(end)+1:end,6:8) = nan;

end


%% set colors
color1 = [255  69   0;... % orange
           67 110 238;... % light blue
          102 205   0;... % light green 
          238 130 238;... % light violet
          255 193  37;... % golden
          128 128 128     % gray
          ]/255;
     
color2 = [220  20  60;... % red
            0   0 139;... % blue
            0 139   0;... % green 
            85 26 139;... % violet
            238 118 0;... % dark orange
            0   0   0     % black
           ]/255;

     
if length(important_centroid) > length(color1)
    color1 = repmat([128 128 128]/255,length(important_centroid),1); %grey
    color2 = repmat([  0   0   0]/255,length(important_centroid),1); %black
end



%% figure
%  plot intensity vs frequency

% fig = climada_figuresize(0.6, 1.1);

for imp_i = length(important_centroid):-1:1
    
    %% historical data
    h(1) = plot(1./int{imp_i}(:,3),  int{imp_i}(:,1), '.' , 'markersize',8,'color',  color1(imp_i,:));
    hold on
    h(2) = plot(1./int{imp_i}(:,3),  int{imp_i}(:,2), '--', 'markersize',3,'color',  color1(imp_i,:));
    h(3) = plot(1./int{imp_i}(:,8),  int{imp_i}(:,6), 'o:', 'markersize',3,'color',  color1(imp_i,:));
    
    %% probabilistical data
    h(4)     = plot(1./int{imp_i}(:,5),  int{imp_i}(:,4), '.' , 'markersize',8,'color',  color2(imp_i,:));
    g(imp_i) = plot(1./int{imp_i}(:,8),  int{imp_i}(:,7), 'o-', 'markersize',3,'color',  color2(imp_i,:));
end

if check_xlog
    set(gca,'XScale','log');
end
set(gca,'XGrid','on')
set(gca,'YGrid','on')
xlabel('Return period (years)')
switch hazard.peril_ID
    case 'TC'
        ylabel('Wind speed [m s^{-1}]')
    case 'TS'
        ylabel('Surge height [m]')
    case 'TR' | 'MA'
        ylabel('Rainfall [mm]')
    otherwise
        if isfield(hazard,'units')
            ylabel(sprintf('Hazard intensity [%s]',hazard.units))
        else ylabel('Hazard intensity'); end
end
title(titlestr)

%legend
legendstr = {'Data',...
            'Data smoothed'  ,...
            ['Gumbel fit hist., \mu = ',num2str(-p(1,1) ,'%10.2f'),', \sigma = ',(num2str(-p(2,1) ,'%10.2f'))],...
            'Prob.'};        
for imp_i = length(important_centroid):-1:1
    legendstr{4+imp_i} = ['Gumbel fit prob., \mu = ',num2str(-p_(1,imp_i),'%10.2f'),', \sigma = ',(num2str(-p_(2,imp_i),'%10.2f')) ' ' name{imp_i}];
end
legend([h g],legendstr,'location','se')
 
if isempty(check_printplot)
    choice = questdlg('print?','print');
    switch choice
        case 'Yes'
            check_printplot = 1;
        case 'No'
            check_printplot = 0;
        case 'Cancel'
            return
    end
end

if check_printplot %(>=1)
    foldername  = [filesep 'results' filesep 'intensity_return_' [name{:}] '.pdf'];
    print(gcf,'-dpdf',[climada_global.data_dir foldername])
    close
    fprintf('saved 1 FIGURE in folder %s \n', foldername);
end
 
end





%%









