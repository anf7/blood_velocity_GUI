function varargout = mapping_GUI(varargin)
% MAPPING_GUI MATLAB code for mapping_GUI.fig
%      MAPPING_GUI, by itself, creates a new MAPPING_GUI or raises the existing
%      singleton*.
%
%      H = MAPPING_GUI returns the handle to a new MAPPING_GUI or the handle to
%      the existing singleton*.
%
%      MAPPING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPPING_GUI.M with the given input arguments.
%
%      MAPPING_GUI('Property','Value',...) creates a new MAPPING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mapping_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mapping_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mapping_GUI

% Last Modified by GUIDE v2.5 30-Jul-2013 21:14:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mapping_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mapping_GUI_OutputFcn, ...
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


% --- Executes just before mapping_GUI is made visible.
function mapping_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mapping_GUI (see VARARGIN)

handles.homedir = pwd;
handles.filetree = [];
handles.downsample = false;
handles.premask = false;
handles.registervid = false;
handles.dynamic_processing = false;
handles.subpixel = false;
handles.createvid = false;
handles.dispaxes = imagesc(ones(500)); axis off; colormap(gray)
handles.framerate = 24;
handles.pixpermm = 500;
handles.fshift = 5;
handles.seqlength = 64;
handles.seqspacing = 64;
handles.first = 1;
handles.last = 'end';
handles.indexcell = cell(0);
handles.storeindexvect = [];

% Choose default command line output for mapping_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mapping_GUI wait for user response (see UIRESUME)
% uiwait(handles.window);


% --- Outputs from this function are returned to the command line.
function varargout = mapping_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fr_rate_Callback(hObject, eventdata, handles)
% hObject    handle to fr_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.framerate = str2double(get(hObject,'String'));
if isnan(handles.framerate) || handles.framerate <= 0
    set(hObject,'String','24');
    handles.framerate = 24;
end

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of fr_rate as text
%        str2double(get(hObject,'String')) returns contents of fr_rate as a double


% --- Executes during object creation, after setting all properties.
function fr_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pix_per_mm_Callback(hObject, eventdata, handles)
% hObject    handle to pix_per_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.pixpermm = str2double(get(hObject,'String'));
if isnan(handles.pixpermm) || handles.pixpermm <= 0
    set(hObject,'String','500');
    handles.pixpermm = 500;
end

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of pix_per_mm as text
%        str2double(get(hObject,'String')) returns contents of pix_per_mm as a double


% --- Executes during object creation, after setting all properties.
function pix_per_mm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pix_per_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fr_size_um_Callback(hObject, eventdata, handles)
% hObject    handle to fr_size_um (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fr_size_um as text
%        str2double(get(hObject,'String')) returns contents of fr_size_um as a double


% --- Executes during object creation, after setting all properties.
function fr_size_um_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_size_um (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.registervid = true;
    enableval = 'on';
else
    handles.registervid = false;
    enableval = 'off';
    set(handles.creat_video,'Value',0);
    set(handles.sub_pix_shift,'Value',0);
    
    handles.registervid = false;
    handles.subpixel = false;
    handles.createvid = false;
end

set(handles.text4,'Enable',enableval);
set(handles.max_frameshift,'Enable',enableval);
set(handles.sub_pix_shift,'Enable',enableval);
set(handles.creat_video,'Enable',enableval);

guidata(hObject, handles);

% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in sub_pix_shift.
function sub_pix_shift_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.subpixel = true;
else
    handles.subpixel = false;
end

guidata(hObject, handles);


% hObject    handle to sub_pix_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sub_pix_shift


% --- Executes on button press in add_file.
function add_file_Callback(hObject, eventdata, handles)

[newfilename,newpathname] = uigetfile('*.tif','Select multi-page TIFF files...','MultiSelect','on');

if ~iscell(newfilename)
    tempfname = newfilename;
    newfilename = cell(1);
    newfilename{1} = tempfname;
end

tiffs = [];
counter = 1;
for n = 1:length(newfilename)
    if length(imfinfo(strcat(newpathname,newfilename{n}))) > 1
         tiffs{counter} = newfilename{n};
         counter = counter + 1;
    else
        warning OFF BACKTRACE
        warning('Cannot load <%s>: Chose only multi-page TIFF files',newfilename{n})
    end
end

if ~isempty(tiffs)
    addfiletree.path = {newpathname};
    addfiletree.tiffs = tiffs;
    addfiletree.level = 0;

    lengthadd = length(addfiletree);
    if lengthadd > 0
        if isempty(handles.filetree)
            handles.filetree = addfiletree;
        else
            match = false;
            x = length(handles.filetree);
            for n = 1:x
                if isequal(handles.filetree(n).path,addfiletree.path)
                    for p = 1:length(addfiletree.tiffs)
                        handles.filetree(n).tiffs(length(handles.filetree(n).tiffs)+1)...
                            = addfiletree.tiffs(p);
                        match = true;
                    end
                end
            end
            if ~match    
                handles.filetree(x+1).path = addfiletree.path;
                handles.filetree(x+1).tiffs = addfiletree.tiffs;
                handles.filetree(x+1).level = addfiletree.level;
            end
        end
    end
    
    handles.filetree = sort_list(handles.filetree);
    
    [filestr, handles.indexcell] = filetree_string(handles.filetree);
    set(handles.file_list,'String',filestr);
end

guidata(hObject, handles);


% hObject    handle to add_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in add_directory.
function add_directory_Callback(hObject, eventdata, handles)

if ispc
    divider = '\';
else
    divider = '/';
end

newdir = uigetdir(handles.homedir,...
    'Select directory...');
if ~isequal(newdir(end),divider)
    newdir(end+1) = divider;
end

addfiletree.path = {newdir};
addfiletree.level = 0;

cd(newdir)
td = dir('*.tif');
cd(handles.homedir);
tiffs = [];
counter = 1;
for n = 1:length(td)
    if length(imfinfo(strcat(newdir,divider,td(n).name))) > 1
         tiffs{counter} = td(n).name;
         counter = counter + 1;
    end
end

if ~isempty(tiffs)
    addfiletree.tiffs = tiffs;

    lengthadd = length(addfiletree);
    if lengthadd > 0
        if isempty(handles.filetree)
            handles.filetree = addfiletree;
        else
            match = false;
            x = length(handles.filetree);
            for n = 1:x
                if isequal(handles.filetree(n).path,addfiletree.path)
                    for p = 1:length(addfiletree.tiffs)
                        handles.filetree(n).tiffs(length(handles.filetree(n).tiffs)+1)...
                            = addfiletree.tiffs(p);
                    end
                    match = true;
                end
            end
            if ~match   
                handles.filetree(x+1).path = addfiletree.path;
                handles.filetree(x+1).tiffs = addfiletree.tiffs;
                handles.filetree(x+1).level = addfiletree.level;
            end
        end
    end
    
    handles.filetree = sort_list(handles.filetree);
    
    [filestr, handles.indexcell] = filetree_string(handles.filetree);
    set(handles.file_list,'String',filestr);
else
     warning OFF BACKTRACE
     warning('There were no multi-page TIFF files in <%s>',newdir)
end

guidata(hObject, handles);

% hObject    handle to add_directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in file_list.
function file_list_Callback(hObject, eventdata, handles)

if ispc
    divider = '\';
else
    divider = '/';
end

indexvect = get(handles.file_list,'Value');
selectsize = length(indexvect);


if selectsize > 0
    if isempty(setdiff(indexvect,handles.storeindexvect)) && ~isempty(setdiff(handles.storeindexvect,indexvect))

        removevect = setdiff(handles.storeindexvect,indexvect);

        for m = 1:length(removevect)
            indexstr = handles.indexcell{removevect(m)};
            for n = 1:length(handles.indexcell)
                if ~isempty(strfind(handles.indexcell{n},indexstr)) && ~sum(removevect == n)
                    removevect(end+1) = n;
                end     
            end
        end

        indexvect = setdiff(indexvect,removevect);
        handles.storeindexvect = indexvect;

    else

        indexvectadd = setdiff(indexvect,handles.storeindexvect);

        for m = 1:length(indexvectadd)
            indexstr = handles.indexcell{indexvectadd(m)};
            for n = 1:length(handles.indexcell)
                if ~isempty(strfind(handles.indexcell{n},indexstr)) && ~sum(indexvect == n)
                    indexvect(end+1) = n;
                end     
            end
        end

        if selectsize > 1
            for m = 1:length(indexvect)
                if sum(handles.storeindexvect == indexvect(m))
                    handles.storeindexvect = [handles.storeindexvect,indexvect];
                    break
                end
                if m == length(indexvect)
                    handles.storeindexvect = indexvect;
                end
            end
        else
            handles.storeindexvect = indexvect;
        end

    end

    handles.storeindexvect = unique(handles.storeindexvect);

    newindexvect = [];
    for m = 1:length(handles.storeindexvect)
        indexstr1 = handles.indexcell{handles.storeindexvect(m)};
        if ~strcmp(indexstr1(end),divider)
            newindexvect(end+1) = handles.storeindexvect(m);
            continue
        end
        for n = 1:length(handles.storeindexvect)
            indexstr2 = handles.indexcell{handles.storeindexvect(n)};
            if strcmp(indexstr2(end),divider)
                continue
            end  
            if ~isempty(strfind(indexstr2,indexstr1))
                newindexvect(end+1) = handles.storeindexvect(m);
                break
            end
        end
    end

    handles.storeindexvect = newindexvect;

    set(handles.file_list,'Value',handles.storeindexvect);

    guidata(hObject, handles);
end
    

% hObject    handle to file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% Hints: contents = cellstr(get(hObject,'String')) returns file_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from file_list


% --- Executes during object creation, after setting all properties.
function file_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_frameshift_Callback(hObject, eventdata, handles)
% hObject    handle to max_frameshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.fshift = round(str2double(get(hObject,'String')));
if isnan(handles.fshift) || handles.fshift <= 0
    set(hObject,'String','5');
    handles.fshift = 5;
end

guidata(hObject, handles);


% Hints: get(hObject,'String') returns contents of max_frameshift as text
%        str2double(get(hObject,'String')) returns contents of max_frameshift as a double


% --- Executes during object creation, after setting all properties.
function max_frameshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_frameshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in creat_video.
function creat_video_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.createvid = true;
else
    handles.createvid = false;
end

guidata(hObject, handles);

% hObject    handle to creat_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of creat_video


% --- Executes on button press in downsample.
function downsample_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.downsample = true;
else
    handles.downsample = false;
end

guidata(hObject, handles);

% hObject    handle to downsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of downsample


% --- Executes on button press in add_dir_plus_sub.
function add_dir_plus_sub_Callback(hObject, eventdata, handles)

addfiletree = find_subdirectories(handles.homedir);
lengthadd = length(addfiletree);

if lengthadd > 0
    if isempty(handles.filetree)
        handles.filetree = addfiletree;
    else
        x = length(addfiletree);
        y = length(handles.filetree);
        for m = 1:x
            match = false;
            for n = 1:y
                if isequal(handles.filetree(n).path,addfiletree(m).path)
                    handles.filetree(n).level = max(addfiletree(m).level,handles.filetree(n).level);
                    for p = 1:length(addfiletree(m).tiffs)
                        handles.filetree(n).tiffs(length(handles.filetree(n).tiffs)+1)...
                            = addfiletree(m).tiffs(p);
                       
                    end
                    match = true;
                end
            end
            if ~match   
                z = length(handles.filetree);
                handles.filetree(z+1).path = addfiletree(m).path;
                handles.filetree(z+1).tiffs = addfiletree(m).tiffs;
                handles.filetree(z+1).level = addfiletree(m).level;
            end
        end
    end
end


handles.filetree = sort_list(handles.filetree);

[filestr, handles.indexcell] = filetree_string(handles.filetree);
set(handles.file_list,'String',filestr);

guidata(hObject, handles);

% hObject    handle to add_dir_plus_sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in process_selected.
function process_selected_Callback(hObject, eventdata, handles)

indexvect = get(handles.file_list,'Value');
indexstr = cell(0);
count = 1;
for m = 1:length(indexvect)
    if length(handles.indexcell{indexvect(m)}) > 4 && strcmp(handles.indexcell{indexvect(m)}(end-3:end),'.tif')
        indexstr{m} = handles.indexcell{indexvect(m)};   
        count = count + 1;
    else
        warning OFF BACKTRACE
        warning('<%s> is a directory.  Ignoring selection...',handles.indexcell{indexvect(m)})
    end
end



if count ~= 1
    instruct.downsample = handles.downsample;
    instruct.registervid = handles.registervid;
    instruct.dynamic_processing = handles.dynamic_processing;
    instruct.subpixel = handles.subpixel;
    instruct.createvid = handles.createvid;
    instruct.framerate = handles.framerate;
    instruct.pixpermm = handles.pixpermm;
    instruct.fshift = handles.fshift;
    instruct.seqlength = handles.seqlength;
    instruct.seqspacing = handles.seqspacing;
    instruct.first = handles.first;
    instruct.last = handles.last;
    instruct.premask = handles.premask;
    
    [savefilename, savepathname] = uiputfile('*.mat', 'Save data as...');
    for n = 1:length(indexstr)
        instruct.filename = indexstr{n};
        [movement_matrix, togetherimg_p] = process_data(instruct);
        
        
        %%..............
    end
else
    warning OFF BACKTRACE
    warning('No multi-page TIFF files have been selected.')

end




% hObject    handle to process_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in process_all.
function process_all_Callback(hObject, eventdata, handles)
% hObject    handle to process_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in dynamic.
function dynamic_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.dynamic_processing = true;
    enableval = 'on';
else
    handles.dynamic_processing = false;
    enableval = 'off';
end

set(handles.text5,'Enable',enableval);
set(handles.text6,'Enable',enableval);
set(handles.dynamic_seq_length,'Enable',enableval);
set(handles.dynamic_seq_spacing,'Enable',enableval);

guidata(hObject, handles);

% hObject    handle to dynamic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dynamic



function dynamic_seq_length_Callback(hObject, eventdata, handles)
% hObject    handle to dynamic_seq_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.seqlength = round(str2double(get(hObject,'String')));
if isnan(handles.seqlength) || handles.seqlength < 16
    set(hObject,'String','64');
    handles.seqlength = 64;
end

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of dynamic_seq_length as text
%        str2double(get(hObject,'String')) returns contents of dynamic_seq_length as a double


% --- Executes during object creation, after setting all properties.
function dynamic_seq_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dynamic_seq_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dynamic_seq_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to dynamic_seq_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.seqspacing = round(str2double(get(hObject,'String')));
if isnan(handles.seqspacing) || handles.seqspacing <= 0
    set(hObject,'String','64');
    handles.seqspacing = 64;
end

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of dynamic_seq_spacing as text
%        str2double(get(hObject,'String')) returns contents of dynamic_seq_spacing as a double


% --- Executes during object creation, after setting all properties.
function dynamic_seq_spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dynamic_seq_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in delete_elements.
function delete_elements_Callback(hObject, eventdata, handles)


indexvect = get(handles.file_list,'Value');
newfiletree = cell(0);
maxlevel = 0;

for m = 1:length(indexvect)
    indexstr = handles.indexcell{indexvect(m)};
    count1 = 1;
    for n = 1:length(handles.filetree)
        if strfind(handles.filetree(n).path{1},indexstr)
            continue
        else
            newfiletree(count1).path = handles.filetree(n).path;
            newfiletree(count1).level = handles.filetree(n).level;
            maxlevel = max([maxlevel,handles.filetree(n).level]);
            newfiletree(count1).tiffs = cell(0); 
            count2 = 1;
            for p = 1:length(handles.filetree(n).tiffs)
                combstr = strcat(handles.filetree(n).path{1},handles.filetree(n).tiffs{p});
                if strfind(combstr,indexstr)
                    continue
                else
                    newfiletree(count1).tiffs(count2) = handles.filetree(n).tiffs(p);
                    count2 = count2 + 1;
                end
            end
            count1 = count1 + 1; 
        end
    end
    handles.filetree = newfiletree;
end
set(handles.file_list,'Value',1);
handles.storeindexvect = 1;
handles.filetree = sort_list(handles.filetree);

for p = 1:maxlevel+1
    newfiletree = handles.filetree;
    handles.filetree = cell(0);
    counter = 1;
    for m = 1:length(newfiletree)
        if isempty(newfiletree(m).tiffs)
            for n = 1:length(newfiletree)
                if m == n
                    continue
                end
                if ~isempty(strfind(newfiletree(n).path{1},newfiletree(m).path{1}))
                    handles.filetree(counter).path = newfiletree(m).path;
                    handles.filetree(counter).level = newfiletree(m).level;
                    handles.filetree(counter).tiffs = newfiletree(m).tiffs; 
                    counter = counter + 1;
                    break
                end
            end
        else
            handles.filetree(counter).path = newfiletree(m).path;
            handles.filetree(counter).level = newfiletree(m).level;
            handles.filetree(counter).tiffs = newfiletree(m).tiffs; 
            counter = counter + 1;
        end
    end
end
   

[filestr, handles.indexcell] = filetree_string(handles.filetree);
set(handles.file_list,'String',filestr);


guidata(hObject, handles);               
                



% hObject    handle to delete_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function first_frm_Callback(hObject, eventdata, handles)

handles.first = round(str2double(get(hObject,'String')));
if isnan(handles.first) || handles.first <= 0
    set(hObject,'String','1');
    handles.first = 1;
end

if ~strcmp(handles.last,'end') && handles.last < handles.first + 16
    handles.last = handles.first + 16;
    set(handles.last_frm,'String',num2str(handles.last));
end    

guidata(hObject, handles);

% hObject    handle to first_frm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of first_frm as text
%        str2double(get(hObject,'String')) returns contents of first_frm as a double


% --- Executes during object creation, after setting all properties.
function first_frm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to first_frm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function last_frm_Callback(hObject, eventdata, handles)
% hObject    handle to last_frm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



if strcmp(lower(get(hObject,'String')),'end')
    handles.last = 'end';
else
    handles.last = round(str2double(get(hObject,'String')));
    if isnan(handles.last) || handles.last < handles.first + 16
        set(hObject,'String','End');
        handles.last = 'end';
    end
end

guidata(hObject, handles);


% Hints: get(hObject,'String') returns contents of last_frm as text
%        str2double(get(hObject,'String')) returns contents of last_frm as a double


% --- Executes during object creation, after setting all properties.
function last_frm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to last_frm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function mode_Callback(hObject, eventdata, handles)
% hObject    handle to mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function analysis_mode_Callback(hObject, eventdata, handles)

set(hObject,'Checked','on')
set(handles.process_mode,'Checked','off')
handles.mode = 'analysis';


guidata(hObject, handles); 

% hObject    handle to analysis_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function process_mode_Callback(hObject, eventdata, handles)

set(hObject,'Checked','on')
set(handles.analysis_mode,'Checked','off')
handles.mode = 'process';

    
guidata(hObject, handles); 

% hObject    handle to process_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in premask.
function premask_Callback(hObject, eventdata, handles)
% hObject    handle to premask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    handles.premask = true;
else
    handles.premask = false;
end

guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of premask
