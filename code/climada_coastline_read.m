
function coastline = climada_coastline_read(coastline_file, check_plot)

%read the coastline file, data source NOAA
% http://www.ngdc.noaa.gov/mgg/coast/
% NAME:
%   climada_coastline_read
% PURPOSE:
%   read the coastline and have coastline as output
%   within: climada_centroids_distance_to_coast
% CALLING SEQUENCE:
%   coastline = climada_coastline_read(coastline_file)
% EXAMPLE:
%   coastline = climada_coastline_read
% INPUTS:
%   none, if coastline_file empty default file is loaded from globalGDP
%   modul data folder
% OPTIONAL INPUT PARAMETERS:
%   check_plot to create figure
% OUTPUTS:
%   coastline
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, 20121203
%-

if ~exist('coastline_file', 'var'), coastline_file = [] ; end
if ~exist('check_plot'    , 'var'), check_plot     = [] ; end

% set modul data directory
modul_data_dir = [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

if isempty(coastline_file)
    coastline_file = [modul_data_dir filesep 'coastline.txt'];
    if ~exist(coastline_file,'file')
        fprintf('Coastline file not availabe in modul data directory. Unable to proceed.\n')
        coastline = [];
        return
    %else 
    end
end

% read the coastline
fid = fopen(coastline_file);
pts = fscanf(fid,'%f %f\n',[2 inf]);
fclose(fid);
coastline.lon = pts(1,:)';
coastline.lat = pts(2,:)';

if check_plot
    climada_figuresize(0.5,0.8)
    plot(coastline.lon,coastline.lat)
end

% % restructure into different closed polygons
% nan_index = find(isnan(coastline.lon));
% for p_i = 1:length(nan_index)-1
%     coastline_polygon(p_i).lon = [coastline.lon(nan_index(p_i)+1 : nan_index(p_i+1)-1); coastline.lon(nan_index(p_i)+1)];
%     coastline_polygon(p_i).lat = [coastline.lat(nan_index(p_i)+1 : nan_index(p_i+1)-1); coastline.lat(nan_index(p_i)+1)];
% end
% 
% % coastline.lon(nan_index) = [];
% % coastline.lat(nan_index) = [];

% figure
% plot(coastline.lon, coastline.lat,'.-')
% hold on
% axis equal
% nan_index = find(isnan(coastline.lon));
% for i = 1:length(nan_index)
%     plot(coastline.lon(nan_index(i)+1:nan_index(i+1)-1),...
%          coastline.lat(nan_index(i)+1:nan_index(i+1)-1),'-rs')
%     pause
% end






