function varargout = climada_track_map(varargin)
% CLIMADA_TRACK_MAP M-file for climada_track_map.fig
%      CLIMADA_TRACK_MAP, by itself, creates a new CLIMADA_TRACK_MAP or raises the existing
%      singleton*.
%
%      H = CLIMADA_TRACK_MAP returns the handle to a new CLIMADA_TRACK_MAP or the handle to
%      the existing singleton*.
%
%      CLIMADA_TRACK_MAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLIMADA_TRACK_MAP.M with the given input arguments.
%
%      CLIMADA_TRACK_MAP('Property','Value',...) creates a new CLIMADA_TRACK_MAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before climada_track_map_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to climada_track_map_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help climada_track_map

% Last Modified by GUIDE v2.5 08-Feb-2013 18:15:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @climada_track_map_OpeningFcn, ...
                   'gui_OutputFcn',  @climada_track_map_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before climada_track_map is made visible.
function climada_track_map_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to climada_track_map (see VARARGIN)

% Choose default command line output for climada_track_map
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using climada_track_map.
if strcmp(get(hObject,'Visible'),'off')
    %plot(rand(5));
end
global climada_global
global tc_track; global hazard; global hazard_2030; 
global centroids; global ELS; global entity; 
global asset_handles; global cmap; global cbar_assets; global miv; global mav

region_str = {'USA' 'China' 'Japan' 'Australia'};
% region_str = {'USA' 'China' 'Japan'};
set(handles.popupmenu6, 'String', region_str);
% popupmenu6_Callback(hObject, eventdata, handles)

set(handles.checkbox3,'String',region_str{1},'Value',0);
set(handles.checkbox4,'String',region_str{2});
set(handles.checkbox5,'String',region_str{3});
set(handles.checkbox6,'String',region_str{4});

% entity
fprintf('Load assets: %s... ', sprintf(' %s',region_str{:}))
modul_data_dir = [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
folder_name    = [modul_data_dir filesep 'track_map' filesep];
clear assets
for r_i = 1:length(region_str)
    entity_file    = [folder_name region_str{r_i} '_entity_2012'];
    load(entity_file)
    assets(r_i) = entity.assets;
end
fprintf('done.\n')

% % display world in ax 1
axes(handles.axes1);
cla;
% % climada_plot_world_borders(0.7,'United States (USA)');
% % climada_plot_world_borders(0.7,'China');
climada_plot_world_borders(0.7);
hold on
xlabel('Longitude')
ylabel('Latitude')
set(gca,'layer','top')

% plot the assets
markersize    = 2;
asset_handles = {};
% miv = min(arrayfun(@(x) (min(x.Value)),assets));
% mav = max(arrayfun(@(x) (max(x.Value)),assets));
% miv = []; mav = [];
beginColor  = [232 232 232 ]/255;
middleColor = [105 105 105 ]/255;
cmap1 = makeColorMap(beginColor, middleColor, 4);
cmap2 = makeColorMap([255 236 139]/255, [255 97 3 ]/255, 6); %[255 153 18]/255 yellow
cmap3 = makeColorMap([205 150 205 ]/255, [93 71 139 ]/255, 2); %[255 153 18]/255 yellow
cmap  = [cmap1; cmap2; cmap3];

for r_i = 1:length(region_str)
    [cbar asset_handles{r_i}]= plotclr(assets(r_i).Longitude, assets(r_i).Latitude, assets(r_i).Value, 's',markersize, 0,0,[],cmap,1,0); 
    miv(r_i) = min(assets(r_i).Value);
    mav(r_i) = max(assets(r_i).Value);
end
set([asset_handles{:}],'markersize',3)
set([asset_handles{:}],'HandleVisibility','off')

axes(handles.axes5);
axis off
colormap(cmap)
% freezeColors(handles.axes5) %freezeColors(cbar_assets)
cbar_assets = colorbar('location','north'); %southoutside
caxis([0 1])
set(cbar_assets,'xlim',[0 1],'xtick',[])
cbar_assets = cbfreeze(cbar_assets);

axes(handles.axes4);
modul_data_dir = [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
folder_name    = [modul_data_dir filesep 'track_map' filesep];
bmp_file = [folder_name 'saffir_simpson_scale.bmp'];
% I = imread(bmp_file);
% imshow(imresize(I, 3))
% imshow(bmp_file)

checkbox3_Callback(hObject, eventdata, handles)
checkbox4_Callback(hObject, eventdata, handles)
checkbox5_Callback(hObject, eventdata, handles)
checkbox6_Callback(hObject, eventdata, handles)

% % axis equal
% set(gca,'DataAspectRatio',[1 1 1]);
% hb = get(gca,'PlotBoxAspectRatio');
% hb = hb(1)/hb(2);

axes(handles.axes2);
% ylabel('Wind speed (kn)')
% xlabel('Time (h)')
% 
% tc_years = sort(unique([tc_track(:).yyyy]),'descend');
% set(handles.popupmenu3, 'String', tc_years);
% 
% tc_cat  = sort(unique([tc_track(:).category]),'descend');
% set(handles.popupmenu5, 'String', tc_cat);
% 
% popupmenu3_Callback(hObject, eventdata, handles)
% popupmenu4_Callback(hObject, eventdata, handles)
% 
% if size(hazard.arr,1) == length(tc_track)
%     tc_maxwind    = full(max(hazard.arr'));
%     [tc_maxwind track_nos] = sort(tc_maxwind,'descend');
%     track_nos = track_nos(1:500);
%     % check for historical tracks only
%     hist_only = get(handles.checkbox2, 'Value');
%     if hist_only == 1
%         ori_flag      = [tc_track.orig_event_flag];
%         track_nos     = track_nos(logical(ori_flag(track_nos)));
%     end
%     listbox5_str  = {length(track_nos)};
%     for t_i = 1:length(track_nos)
%         listbox5_str{t_i} = sprintf('%d, \t %s - %s, \t\t %d \t\t %2.1f m/s, %s', track_nos(t_i), ...
%             datestr(tc_track(track_nos(t_i)).nodetime_mat(1),'dd/mm'), ...
%             datestr(tc_track(track_nos(t_i)).nodetime_mat(end),'dd/mm/yyyy'), ...
%             tc_track(track_nos(t_i)).category, tc_maxwind(t_i), strrep(tc_track(track_nos(t_i)).name,' ',''));
%     end
%     set(handles.listbox5, 'String', listbox5_str);
% end


% UIWAIT makes climada_track_map wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = climada_track_map_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tc_track
global asset_handles; global cbar_assets
% set([asset_handles{:}],'HandleVisibility','off')

popup_sel_index = 1;
% popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        
        no_hist      = sum([tc_track.orig_event_flag]);
        no_generated = length(tc_track)/no_hist;
        ens_size     = no_generated-1;        
        track_no     = str2num(get(handles.edit1, 'String'));
        if track_no>length(tc_track)
            return
        end
        axes(handles.axes1);
        cla
        
        if get(handles.checkbox9,'Value')
            x_range = [min(tc_track(track_no).lon)-5 max(tc_track(track_no).lon)+5];
            y_range = [min(tc_track(track_no).lat)-5 max(tc_track(track_no).lat)+5];

            hb = get(gca,'PlotBoxAspectRatio');
            hb = hb(1)/hb(2);

            if hb/(diff(x_range)/diff(y_range))<1
                dif     = ( diff(x_range)/hb-diff(y_range) )/2;
                y_range = [y_range(1)-dif y_range(2)+dif];
                set(gca,'xlim',x_range,'ylim',y_range)
            else
                dif     = ( diff(y_range)*hb-diff(x_range) )/2;
                x_range = [x_range(1)-dif x_range(2)+dif]; 
                set(gca,'xlim',x_range,'ylim',y_range)
            end
            climada_plot_world_borders(0.7,[], [], 1)

            tc_track_.lon = tc_track(track_no).lon(1:6:end);
            tc_track_.lat = tc_track(track_no).lat(1:6:end);
            tc_track_.MaxSustainedWind = tc_track(track_no).MaxSustainedWind(1:6:end);
            climada_plot_tc_track_stormcategory(tc_track_,8,[]);
            plot(tc_track(track_no).lon(1:6:end),tc_track(track_no).lat(1:6:end),'ok','markersize',3,'linewidth',0.7);
            % g = plot(tc_track(track_req).lon,tc_track(track_req).lat,'ok','markersize',3,'linewidth',0.7);
            title(['Track ' int2str(track_no) ', single track only'])
            
            axes(handles.axes2);
            cla
            track_no = climada_plot_probabilistic_wind_speed_decay(tc_track, track_no, 1);
        else
            climada_plot_probabilistic_wind_speed_map_gui(tc_track, track_no);
            %[track_no h g] = climada_plot_probabilistic_wind_speed_map_gui(tc_track, track_no);

            axes(handles.axes2);
            cla
            track_no = climada_plot_probabilistic_wind_speed_decay(tc_track, track_no);
        end
        
        %set(handles.edit1, 'String',int2str(track_no+ens_size+1))
        
        checkbox3_Callback(hObject, eventdata, handles)
        checkbox4_Callback(hObject, eventdata, handles)
        checkbox5_Callback(hObject, eventdata, handles)
        checkbox6_Callback(hObject, eventdata, handles)
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
if get(handles.radiobutton1, 'Value')
    pushbutton1_Callback(hObject, eventdata, handles)
elseif get(handles.radiobutton2, 'Value')
    pushbutton4_Callback(hObject, eventdata, handles)
elseif get(handles.radiobutton3, 'Value')
    pushbutton6_Callback(hObject, eventdata, handles)    
end



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tc_track
global hazard; global hazard_2030
track_no = str2double(get(handles.edit1, 'String'));
if track_no>length(tc_track)
    return
end

clim = get(handles.checkbox10,'Value');
if clim == 1 && ~isempty(hazard_2030)
    hazard_.arr = hazard_2030.arr(track_no,:);
    hazard_.lon = hazard_2030.lon;
    hazard_.lat = hazard_2030.lat;
    nametag     = sprintf('%d (cc 2030)',track_no);
else
    hazard_.arr = hazard.arr(track_no,:);
    hazard_.lon = hazard.lon;
    hazard_.lat = hazard.lat;
    nametag     = int2str(track_no);
end

axes(handles.axes1);
% set([asset_handles{:}],'HandleVisibility','off')
cla

climada_plot_footprint(hazard_, tc_track(track_no), nametag);

% climada_plot_windfield(hazard, tc_track, track_no);
% % [contr t_handle] = climada_plot_windfield(hazard, tc_track, track_no);

axes(handles.axes2);
cla

checkbox3_Callback(hObject, eventdata, handles)
checkbox4_Callback(hObject, eventdata, handles)
checkbox5_Callback(hObject, eventdata, handles)
checkbox6_Callback(hObject, eventdata, handles)


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
global tc_track
year_all  = str2num(get(handles.popupmenu3, 'String'));
year_sel  = year_all(get(handles.popupmenu3, 'Value'));
tc_years  = arrayfun(@(x) (min(x.yyyy)), tc_track);
tc_index  = tc_years == year_sel;
% tc_track_ = tc_track(tc_index);
tc_month  = arrayfun(@(x) (min(x.mm)), tc_track(tc_index));
set(handles.popupmenu4, 'String', unique(tc_month));

track_nos = find(tc_index);
% check for historical tracks only
ori_flag      = [tc_track.orig_event_flag];
track_nos     = track_nos(logical(ori_flag(track_nos)));
cat_sel       = [tc_track(track_nos).category];
tc_startdates = arrayfun(@(x) (min(x.nodetime_mat)), tc_track(track_nos));
tc_enddates   = arrayfun(@(x) (max(x.nodetime_mat)), tc_track(track_nos));
listbox1_str  = {};
for t_i = 1:length(tc_startdates)
    listbox1_str{t_i} = sprintf('%d, \t %s - %s, \t\t %d', track_nos(t_i), datestr(tc_startdates(t_i),'dd/mm'), datestr(tc_enddates(t_i),'dd/mm'), cat_sel(t_i));
end
set(handles.listbox1, 'String', listbox1_str);



% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
global tc_track
year_all  = str2num(get(handles.popupmenu3, 'String'));
year_sel  = year_all(get(handles.popupmenu3, 'Value'));
month_all = str2num(get(handles.popupmenu4, 'String'));
month_sel = month_all(get(handles.popupmenu4, 'Value'));

tc_years  = arrayfun(@(x) (min(x.yyyy)), tc_track);
tc_index  = tc_years == year_sel;
tc_track_ = tc_track(tc_index);
tc_month  = arrayfun(@(x) (min(x.mm)), tc_track_);
tc_index2 = tc_month == month_sel;
t         = find(tc_index);
track_nos = t(tc_index2);
% check for historical tracks only
ori_flag      = [tc_track.orig_event_flag];
track_nos     = track_nos(logical(ori_flag(track_nos)));
cat_sel       = [tc_track(track_nos).category];
tc_startdates = arrayfun(@(x) (min(x.nodetime_mat)), tc_track(track_nos));
tc_enddates   = arrayfun(@(x) (max(x.nodetime_mat)), tc_track(track_nos));
listbox1_str  = {};
for t_i = 1:length(tc_startdates)
    listbox1_str{t_i} = sprintf('%d, \t %s - %s, \t\t %d', track_nos(t_i), datestr(tc_startdates(t_i),'dd/mm'), datestr(tc_enddates(t_i),'dd/mm'), cat_sel(t_i));
end
set(handles.listbox1, 'String', listbox1_str);




% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
listbox1_str  = get(handles.listbox1, 'String');
listbox1_str  = listbox1_str(get(handles.listbox1, 'Value'));
track_no      = textscan(listbox1_str{:}, '%d %s %d', 'delimiter', ',');
track_no      = track_no{1};
set(handles.edit1, 'String',track_no)
if get(handles.radiobutton1, 'Value')
    pushbutton1_Callback(hObject, eventdata, handles)
elseif get(handles.radiobutton2, 'Value')
    pushbutton4_Callback(hObject, eventdata, handles)
elseif get(handles.radiobutton3, 'Value')
    pushbutton6_Callback(hObject, eventdata, handles)    
end



% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on popupmenu3 and none of its controls.
function popupmenu3_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4
listbox4_str  = get(handles.listbox4, 'String');
listbox4_str  = listbox4_str(get(handles.listbox4, 'Value'));
track_no      = textscan(listbox4_str{:}, '%d %s', 'delimiter', ',');
track_no      = track_no{1};
set(handles.edit1, 'String',track_no)
if get(handles.radiobutton1, 'Value')
    pushbutton1_Callback(hObject, eventdata, handles)
elseif get(handles.radiobutton2, 'Value')
    pushbutton4_Callback(hObject, eventdata, handles)
elseif get(handles.radiobutton3, 'Value')
    pushbutton6_Callback(hObject, eventdata, handles) 
end


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
global tc_track
cat_all   = str2num(get(handles.popupmenu5, 'String'));
cat_sel   = cat_all(get(handles.popupmenu5, 'Value'));
tc_cat    = [tc_track(:).category];
track_nos = find(tc_cat == cat_sel);
% % check for historical tracks only
hist_only = get(handles.checkbox1, 'Value');
if hist_only == 1
    ori_flag      = [tc_track.orig_event_flag];
	track_nos     = track_nos(logical(ori_flag(track_nos)));
end
tc_startdates = arrayfun(@(x) (min(x.nodetime_mat)), tc_track(track_nos));
tc_enddates   = arrayfun(@(x) (max(x.nodetime_mat)), tc_track(track_nos));
listbox4_str  = {};
for t_i = 1:length(tc_startdates)
    listbox4_str{t_i} = sprintf('%d, \t %s - %s', track_nos(t_i), datestr(tc_startdates(t_i),'dd/mm'), datestr(tc_enddates(t_i),'dd/mm/yyyy'));
end
set(handles.listbox4, 'String', listbox4_str);




% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
global tc_track
cat_all   = str2num(get(handles.popupmenu5, 'String'));
cat_sel   = cat_all(get(handles.popupmenu5, 'Value'));
tc_cat    = [tc_track(:).category];
track_nos = find(tc_cat == cat_sel);
% % check for historical tracks only
hist_only = get(handles.checkbox1, 'Value');
if hist_only == 1
    ori_flag      = [tc_track.orig_event_flag];
	track_nos     = track_nos(logical(ori_flag(track_nos)));
end
tc_startdates = arrayfun(@(x) (min(x.nodetime_mat)), tc_track(track_nos));
tc_enddates   = arrayfun(@(x) (max(x.nodetime_mat)), tc_track(track_nos));
listbox4_str  = {};
for t_i = 1:length(tc_startdates)
    listbox4_str{t_i} = sprintf('%d, \t %s - %s', track_nos(t_i), datestr(tc_startdates(t_i),'dd/mm'), datestr(tc_enddates(t_i),'dd/mm/yyyy'));
end
set(handles.listbox4, 'String', listbox4_str);


% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5
listbox5_str  = get(handles.listbox5, 'String');
listbox5_str  = listbox5_str(get(handles.listbox5, 'Value'));
track_no      = textscan(listbox5_str{:}, '%d %s %s', 'delimiter', ',');
track_no      = track_no{1};
set(handles.edit1, 'String',track_no)
if get(handles.radiobutton1, 'Value')
    pushbutton1_Callback(hObject, eventdata, handles)
elseif get(handles.radiobutton2, 'Value')
    pushbutton4_Callback(hObject, eventdata, handles)
elseif get(handles.radiobutton3, 'Value')
    pushbutton6_Callback(hObject, eventdata, handles) 
end


% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
global tc_track
global hazard
tc_maxwind    = full(max(hazard.arr'));
[tc_maxwind track_nos] = sort(tc_maxwind,'descend');
track_nos = track_nos(1:448);
% check for historical tracks only
hist_only = get(handles.checkbox2, 'Value');
if hist_only == 1
    ori_flag      = [tc_track.orig_event_flag];
	track_nos     = track_nos(logical(ori_flag(track_nos)));
end
listbox5_str  = {length(track_nos)};

for t_i = 1:length(track_nos)
    listbox5_str{t_i} = sprintf('%d, \t %s - %s, \t\t %d \t\t %2.1f m/s, %s', track_nos(t_i), ...
        datestr(tc_track(track_nos(t_i)).nodetime_mat(1),'dd/mm'), ...
        datestr(tc_track(track_nos(t_i)).nodetime_mat(end),'dd/mm/yyyy'), ...
        tc_track(track_nos(t_i)).category, tc_maxwind(t_i), strrep(tc_track(track_nos(t_i)).name,' ',''));
end
set(handles.listbox5, 'String', listbox5_str);


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1'
        % Code for when radiobutton1 is selected.
        pushbutton1_Callback(hObject, eventdata, handles)
    case 'radiobutton2'
        % Code for when radiobutton2 is selected.
        pushbutton4_Callback(hObject, eventdata, handles)
    case 'radiobutton3'
        % Code for when radiobutton2 is selected.
        pushbutton6_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tc_track; global ELS; global centroids

track_no = str2num(get(handles.edit1, 'String'));
if track_no>length(tc_track)
    return
end

clim = get(handles.checkbox11,'Value');
if clim == 1
    event_loss = ELS(3).loss_per_cu(track_no,:);
    nametag     = sprintf('%d (cc 2030)',track_no);
else
    event_loss = ELS(1).loss_per_cu(track_no,:);
    nametag     = int2str(track_no);
end

axes(handles.axes1);
cla
climada_plot_lossfootprint(event_loss, centroids, tc_track(track_no), nametag)

axes(handles.axes2);
cla
cla(gca, 'reset');
climada_plot_loss_hist(event_loss, nametag)

checkbox3_Callback(hObject, eventdata, handles)
checkbox4_Callback(hObject, eventdata, handles)
checkbox5_Callback(hObject, eventdata, handles)
checkbox6_Callback(hObject, eventdata, handles)




% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6

global climada_global
global hazard; global hazard_2030; global tc_track; global ELS; global centroids
region_all  = get(handles.popupmenu6, 'String');
region_sel  = region_all{get(handles.popupmenu6, 'Value')};

% load tracks, hazard, centroids and ELS
modul_data_dir   = [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
folder_name      = [modul_data_dir filesep 'track_map' filesep];
tc_track_file    = [folder_name region_sel '_tc_tracks.mat'];
hazard_file      = [folder_name region_sel '_hazard_2012.mat'];
hazard_clim_file = [folder_name region_sel '_hazard_2030.mat'];
ELS_file         = [folder_name region_sel '_ELS_2012_2030eco_2030clim.mat'];
centroids_file   = [folder_name region_sel '_centroids.mat'];
if strcmp(region_sel,'Japan')
    region_sel    = 'China';
    tc_track_file = [folder_name region_sel '_tc_tracks.mat'];
end
fprintf('Load tc tracks, hazard and ELS: %s...', region_sel)
tc_track = []; hazard = []; hazard_2030 = [];
load(tc_track_file)
load([folder_name 'Mozambique' '_tc_tracks.mat'])
load(hazard_file); 
% hazard = hazard_wkn;
vars = whos('-file', hazard_file);
if ~strcmp(vars.name,'hazard')
        hazard = eval(vars.name);
        clear (vars.name)
end
if exist(hazard_clim_file,'file') == 2
    load(hazard_clim_file)
    vars = whos('-file', hazard_clim_file);
    if ~strcmp(vars.name,'hazard_2030')
        hazard_2030 = eval(vars.name);
        clear (vars.name)
    end
end
load(ELS_file)
vars = whos('-file', ELS_file);
if ~strcmp(vars.name,'ELS')
        ELS = eval(vars.name);
        clear (vars.name)
end
load(centroids_file)
fprintf(' done.\n')

tc_years = sort(unique([tc_track(:).yyyy]),'descend');
set(handles.popupmenu3, 'String', tc_years);

tc_cat  = sort(unique([tc_track(:).category]),'descend');
set(handles.popupmenu5, 'String', tc_cat);

popupmenu3_Callback(hObject, eventdata, handles)
popupmenu4_Callback(hObject, eventdata, handles)
popupmenu5_Callback(hObject, eventdata, handles)
checkbox2_Callback(hObject, eventdata, handles)

% tc_maxwind    = full(max(hazard.arr'));
% [tc_maxwind track_nos] = sort(tc_maxwind,'descend');
% track_nos = track_nos(1:500);
% % check for historical tracks only
% hist_only = get(handles.checkbox2, 'Value');
% if hist_only == 1
%     ori_flag      = [tc_track.orig_event_flag];
% 	track_nos     = track_nos(logical(ori_flag(track_nos)));
% end
% listbox5_str  = {length(track_nos)};
% for t_i = 1:length(track_nos)
%     listbox5_str{t_i} = sprintf('%d, \t %s - %s, \t\t %d \t\t %2.1f m/s', track_nos(t_i), ...
%         datestr(tc_track(track_nos(t_i)).nodetime_mat(1),'dd/mm'), ...
%         datestr(tc_track(track_nos(t_i)).nodetime_mat(end),'dd/mm/yyyy'), ...
%         tc_track(track_nos(t_i)).category, tc_maxwind(t_i) );
% end
% set(handles.listbox5, 'String', listbox5_str);




% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
global asset_handles; global cbar_assets; global miv; global mav
if get(handles.checkbox3, 'Value')
    set(asset_handles{1},'handlevisibility','on','visible','on')
    set(cbar_assets,'xlim',[0 12],'xtick',[6 11],'xticklabel',{num2str(mav(1)*0.5,'%2.1e') num2str(mav(1),'%2.1e')})
else
    set(asset_handles{1},'handlevisibility','on','visible','off')
end
set(asset_handles{1},'handlevisibility','off')


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox4
global asset_handles; global cbar_assets; global miv; global mav
if get(handles.checkbox4, 'Value')
    set(asset_handles{2},'visible','on')
    set(cbar_assets,'xlim',[0 12],'xtick',[6 11],'xticklabel',{num2str(mav(2)*0.5,'%2.1e') num2str(mav(2),'%2.1e')})
else
    set(asset_handles{2},'visible','off')
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox5
global asset_handles; global cbar_assets; global miv; global mav
if get(handles.checkbox5, 'Value')
    set(asset_handles{3},'visible','on')
    set(cbar_assets,'xlim',[0 12],'xtick',[6 11],'xticklabel',{num2str(mav(3)*0.5,'%2.1e') num2str(mav(3),'%2.1e')})
else
    set(asset_handles{3},'visible','off')
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
global asset_handles; global cbar_assets; global miv; global mav
if get(handles.checkbox6, 'Value')
    set(asset_handles{4},'visible','on')
    set(cbar_assets,'xlim',[0 12],'xtick',[6 11],'xticklabel',{num2str(mav(3)*0.5,'%2.1e') num2str(mav(3),'%2.1e')})
else
    set(asset_handles{4},'visible','off')
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
track_no = str2num(get(handles.edit1, 'String'));
set(handles.edit1, 'String',int2str(track_no+1));
%update graphics
edit1_Callback(hObject, eventdata, handles)



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
track_no = str2num(get(handles.edit1, 'String'));
set(handles.edit1, 'String',int2str(track_no-1));
%update graphics
edit1_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listbox5_val = get(handles.listbox5, 'Value');
max_val      = size(get(handles.listbox5, 'string'),1);
if listbox5_val<max_val
    set(handles.listbox5, 'Value',listbox5_val+1);
    %update graphics
    listbox5_Callback(hObject, eventdata, handles)
end



% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listbox5_val = get(handles.listbox5, 'Value');
if listbox5_val>1
    set(handles.listbox5, 'Value',listbox5_val-1);
    %update graphics
    listbox5_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listbox4_val = get(handles.listbox4, 'Value');
max_val      = size(get(handles.listbox4, 'string'),1);
if listbox4_val<max_val
    set(handles.listbox4, 'Value',listbox4_val+1);
    %update graphics
    listbox4_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listbox4_val = get(handles.listbox4, 'Value');
if listbox4_val>1
    set(handles.listbox4, 'Value',listbox4_val-1);
    %update graphics
    listbox4_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox9
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4
axes(hObject)
modul_data_dir = [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
folder_name    = [modul_data_dir filesep 'track_map' filesep];
bmp_file = [folder_name 'saffir_simpson_scale.bmp'];
% I = imread(bmp_file);
% imshow(imresize(I, 3))
imshow(bmp_file)


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox10
pushbutton4_Callback(hObject, eventdata, handles)



% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11
pushbutton6_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
