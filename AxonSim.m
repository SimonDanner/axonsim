function varargout = AxonSim(varargin)
% AxonSim M-file for AxonSim.fig
%      AxonSim, by itself, creates a new AxonSim or raises the existing
%      singleton*.
%
%      H = AxonSim returns the handle to a new AxonSim or the handle to
%      the existing singleton*.
%
%      AxonSim('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AxonSim.M with the given input arguments.
%
%      AxonSim('Property','Value',...) creates a new AxonSim or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AxonSim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AxonSim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AxonSim

% Last Modified by GUIDE v2.5 17-Apr-2011 13:48:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AxonSim_OpeningFcn, ...
                   'gui_OutputFcn',  @AxonSim_OutputFcn, ...
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


% --- Executes just before AxonSim is made visible.
function AxonSim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AxonSim (see VARARGIN)

% Choose default command line output for AxonSim
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AxonSim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AxonSim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_electrode.
function popupmenu_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_electrode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_electrode
val = get(hObject,'Value');

if val == 4
    set(handles.text4,'String','Stimulation strength [pA]:');
    set(handles.text11,'String','Accuracy [pA]:');
    set(handles.text5,'String','Location of electrode (Node):');
    set(handles.edit_strength,'String','600');
    set(handles.edit_e_pos,'String','15');
else
    set(handles.text11,'String','Accuracy [mA]:');
    set(handles.text4,'String','Stimulation strength [mA]:');
    set(handles.text5,'String','Electrode position [mm] [x y]:');
    set(handles.edit_strength,'String','-.2');
    set(handles.edit_e_pos,'String','[0 1]');
end


% --- Executes during object creation, after setting all properties.
function popupmenu_electrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
set(hObject,'String',{'single','double 1/-1','triple -0.5/1/-0.5','intracellular'});
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dia_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dia as text
%        str2double(get(hObject,'String')) returns contents of edit_dia as a double


% --- Executes during object creation, after setting all properties.
function edit_dia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_strength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_strength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_strength as text
%        str2double(get(hObject,'String')) returns contents of edit_strength as a double


% --- Executes during object creation, after setting all properties.
function edit_strength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_strength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_e_pos_Callback(hObject, eventdata, handles)
% hObject    handle to edit_e_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_e_pos as text
%        str2double(get(hObject,'String')) returns contents of edit_e_pos as a double


% --- Executes during object creation, after setting all properties.
function edit_e_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_e_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dur_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dur as text
%        str2double(get(hObject,'String')) returns contents of edit_dur as a double


% --- Executes during object creation, after setting all properties.
function edit_dur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_frq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frq as text
%        str2double(get(hObject,'String')) returns contents of edit_frq as a double


% --- Executes during object creation, after setting all properties.
function edit_frq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sim_dur_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sim_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sim_dur as text
%        str2double(get(hObject,'String')) returns contents of edit_sim_dur as a double


% --- Executes during object creation, after setting all properties.
function edit_sim_dur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sim_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
    
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA


e_pos = str2num(get(handles.edit_e_pos,'String'));
c = str2double(get(handles.edit_curve,'String'));
phi = str2double(get(handles.edit_res,'String'));
dia = str2double(get(handles.edit_dia,'String'));
sim_dur = str2double(get(handles.edit_sim_dur,'String'));
pulse_dur  = str2double(get(handles.edit_dur,'String'));
frq  = str2double(get(handles.edit_frq,'String'));
e_sep  = str2double(get(handles.edit_e_sep,'String'));
e_type = get(handles.popupmenu_electrode,'Value');
accuracy = str2double(get(handles.edit_accuracy,'String'));
length = str2double(get(handles.edit_length,'String'));


e_sep = e_sep*1000;
calc_thr = get(handles.check_thr,'Value');
fun_type = get(handles.popupmenu_stim,'Value');
model_nr = get(handles.popupmenu_model,'Value');
custom_fun = get(handles.edit_stimfun,'String');

e_pos = e_pos * 1000;  %um
phi = phi*100000; % ohm um

I = str2double(get(handles.edit_strength,'String'));
if e_type ~= 4
    I = I /1000; %A
    accuracy = accuracy/1000;
end



[x,y] = curvature(c,length);
l=0;

s = size(x);
for i = 2:s(2)
   l = [l,l(i-1)+sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2)]; 
end


if calc_thr == 1
    I = sign(I) * accuracy;
end


V=[];
switch e_type
    case 1
        V=electrode(phi,I,e_pos(1),e_pos(2),x,y);
        plot(handles.axes1,e_pos(1),e_pos(2),'x'); 
        hold(handles.axes1,'on');
    case 2
        V=electrode(phi,I,e_pos(1)-e_sep/2,e_pos(2),x,y)-electrode(phi,I,e_pos(1)+e_sep/2,e_pos(2),x,y);   
        plot(handles.axes1,e_pos(1)-e_sep/2,e_pos(2),'x'); 
        hold(handles.axes1,'on');
        plot(handles.axes1,e_pos(1)+e_sep/2,e_pos(2),'x'); 
        hold(handles.axes1,'on');
    case 3
        V=electrode(phi,I,e_pos(1),e_pos(2),x,y)-0.5*electrode(phi,I,e_pos(1)-e_sep,e_pos(2),x,y)-0.5*electrode(phi,I,e_pos(1)+e_sep,e_pos(2),x,y);
        plot(handles.axes1,e_pos(1)-e_sep/2,e_pos(2),'x'); 
        hold(handles.axes1,'on');
        plot(handles.axes1,e_pos(1)+e_sep/2,e_pos(2),'x'); 
        hold(handles.axes1,'on');
        plot(handles.axes1,e_pos(1),e_pos(2),'x'); 
        hold(handles.axes1,'on'); 
    case 4
        V=zeros(1,s(2));
end




x_m = x'/1000000;
y_m = y'/1000000;
z_m = zeros(s(2),s(1));
data = [x_m,y_m,z_m,V'];

sI = sign(I);
I = abs(I);
I_orig = I;
V_orig = V;
last_super = I_orig;
resultt=[];
if calc_thr == 1
    step = accuracy/2;
    bigger = 0;
    
    while 1
        if e_type == 4
            [t,Y,N,b_thr] = model(model_nr,sim_dur,data,pulse_dur,fun_type,custom_fun,dia,frq,1,sI*I,e_pos/1000+8);
        else
            [t,Y,N,b_thr] = model(model_nr,sim_dur,data,pulse_dur,fun_type,custom_fun,dia,frq,1);
        end
        
        if b_thr == 1 
            step = step / 2;
            last_super = I;
            I = I - step;
            bigger = 1;
        else
            if bigger == 1
                step = step / 2;
            else
                step = step * 2;
            end
            I = I + step;
        end
        if e_type == 4
            I_pa = I *sI
        else
            I_mA=I*sI*1000
        end
        if step < accuracy
            if e_type~=4
                resultt=last_super*sI*1000;
            else
                resultt=last_super*sI; 
            end
            set(handles.text_thr,'String',resultt);
           break; 
        end
        V = V_orig.*(I/I_orig);
        
        data = [x_m,y_m,z_m,V'];
    end
    V = V_orig.*(last_super/I_orig);
    data = [x_m,y_m,z_m,V'];
end

    
if e_type == 4
    [t,Y,N] = model(model_nr,sim_dur,data,pulse_dur,fun_type,custom_fun,dia,frq,0,sI*I,e_pos/1000+8);
else
    [t,Y,N] = model(model_nr,sim_dur,data,pulse_dur,fun_type,custom_fun,dia,frq,0);
end



af = [];
for i = 2:s(2)-1
    af(i-1) = (V(i)-V(i-1))/(l(i)-l(i-1))^2+(V(i)-V(i+1))/(l(i)-l(i+1))^2;
end

szy = size(Y);

    m = 16;

Yp=zeros(szy(1),N-m);
for i = 1:N-m
   Yp(:,i) = Y(:,i+m/2)-(i-1)*40; 
end


stimf=[];
ts=[];
isi = 1/frq*1000;


for i = 0:1/1000:sim_dur
    stimf=[stimf,sI*stimulation_fun(mod(i,isi),pulse_dur,fun_type,custom_fun)];

    ts=[ts,i];
end

set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:');

plot(handles.axes6,ts,stimf);
xlim(handles.axes6,[0 sim_dur]);
ylim(handles.axes6,[-1.1 1.1]);

plot(handles.axes2,l,V);
xlim(handles.axes2,[0 l(s(2))]);
plot(handles.axes3,l(2:s(2)-1),-af);
xlim(handles.axes3,[0 l(s(2))]);



plot(handles.axes1,x,y)
hold(handles.axes1,'off');
xlim(handles.axes1,[-length/2 length/2]);
if e_type ~= 4
    ylim(handles.axes1,[-length/pi e_pos(2) + 1000]);
end
plot(handles.axes4,t,Y(:,floor(N/2)));
plot(handles.axes5,t,Yp,'-');
xlim(handles.axes4,[0 sim_dur]);
xlim(handles.axes5,[0 sim_dur]);

set(handles.axes5,'YTick',[]);




%ylim(handles.axes1,[-0.1*electrode_distance,1.1*electrode_distance]);



function edit_curve_Callback(hObject, eventdata, handles)
% hObject    handle to edit_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_curve as text
%        str2double(get(hObject,'String')) returns contents of edit_curve as a double


% --- Executes during object creation, after setting all properties.
function edit_curve_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_res_Callback(hObject, eventdata, handles)
% hObject    handle to edit_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_res as text
%        str2double(get(hObject,'String')) returns contents of edit_res as a double


% --- Executes during object creation, after setting all properties.
function edit_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_e_sep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_e_sep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_e_sep as text
%        str2double(get(hObject,'String')) returns contents of edit_e_sep as a double


% --- Executes during object creation, after setting all properties.
function edit_e_sep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_e_sep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_thr.
function check_thr_Callback(hObject, eventdata, handles)
% hObject    handle to check_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_thr


% --- Executes during object creation, after setting all properties.
function text_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit_accuracy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_accuracy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_accuracy as text
%        str2double(get(hObject,'String')) returns contents of edit_accuracy as a double


% --- Executes during object creation, after setting all properties.
function edit_accuracy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_accuracy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stimfun_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stimfun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stimfun as text
%        str2double(get(hObject,'String')) returns contents of edit_stimfun as a double


% --- Executes during object creation, after setting all properties.
function edit_stimfun_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stimfun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_stim.
function popupmenu_stim_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_stim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_stim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_stim


% --- Executes during object creation, after setting all properties.
function popupmenu_stim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_stim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
set(hObject,'String',{'single pulse','double pulse','custom'});
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_length_Callback(hObject, eventdata, handles)
% hObject    handle to edit_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_length as text
%        str2double(get(hObject,'String')) returns contents of edit_length as a double


% --- Executes during object creation, after setting all properties.
function edit_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_model.
function popupmenu_model_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_model


% --- Executes during object creation, after setting all properties.
function popupmenu_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
set(hObject,'String',{'MRG','CRRSS','McNeal'});
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
