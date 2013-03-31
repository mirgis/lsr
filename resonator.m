%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project:   FYO - Laser resonator and stability
% File:      resonator.m
% Date:      30.3.2013
% Author(s): Radek Fér       <xferra00@stud.fit.vutbr.cz>
%            Miroslav Skácel <xskace00@stud.fit.vutbr.cz>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = resonator(varargin)
% RESONATOR M-file for resonator.fig
%      RESONATOR, by itself, creates a new RESONATOR or raises the existing
%      singleton*.
%
%      H = RESONATOR returns the handle to a new RESONATOR or the handle to
%      the existing singleton*.
%
%      RESONATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESONATOR.M with the given input arguments.
%
%      RESONATOR('Property','Value',...) creates a new RESONATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before resonator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to resonator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help resonator

% Last Modified by GUIDE v2.5 31-Mar-2013 23:22:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @resonator_OpeningFcn, ...
                   'gui_OutputFcn',  @resonator_OutputFcn, ...
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


% --- Executes just before resonator is made visible.
function resonator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to resonator (see VARARGIN)

% init
resonator_paint(handles);
stability_paint(handles);
energy_paint(handles);

% Choose default command line output for resonator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes resonator wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Function for painting static data
function resonator_paint(handles)
axes(handles.resonator);
PI=3.1415926535;
L=get(handles.sliderL,'value'); % mirrors distance
r1=str2double(get(handles.editR1,'String')); % M1 radius
r2=str2double(get(handles.editR2,'String')); % M2 radius
t1=-PI:0.001:0; % M1 arc base
t2=0:0.001:PI;  % M2 arc base
x1=sin(t1);
y1=cos(t1);
x2=sin(t2);
y2=cos(t2);
thr1=min(r1/3,0.95); % rendering arc threshold
thr2=min(r2/3,0.95); % rendering arc threshold
xx=x1((r1*y1>-thr1)&(r1*y1<thr1)); % M1 x arc trim
yy=y1((r1*y1>-thr1)&(r1*y1<thr1)); % M1 y arc trim
xxx=x2((r2*y2>-thr2)&(r2*y2<thr2));% M2 x arc trim 
yyy=y2((r2*y2>-thr2)&(r2*y2<thr2));% M2 y arc trim
plot(r1*xx+r1-L/2,r1*yy,'b','LineWidth',2,'Color','g');
hold on;
plot(r1/2-L/2,0,'xb');
line([r1*xx(1)+r1-L/2, r1-L/2],[r1*(yy(1)), 0],'LineStyle',':');
line([r1*xx(end)+r1-L/2, r1-L/2],[r1*(yy(end)), 0],'LineStyle',':');

plot(r2*xxx-r2+L/2,r2*yyy,'b','LineWidth',2);
plot(-r2/2+L/2,0,'xb');
line([r2*xxx(1)-r2+L/2, -r2+L/2],[r2*(yyy(1)), 0],'LineStyle',':');
line([r2*xxx(end)-r2+L/2, -r2+L/2],[r2*(yyy(end)), 0],'LineStyle',':');
line([-100,100],[0,0],'LineStyle',':','Color','k');
axis([min(-3,min(r1*xx+r1-L/2-1)) max(3,max(r2*xxx-r2+L/2+1)) -r1/2 r1/2]);
%axis([-3 3 -3 3]);
hold off;
%set(gca,'xtick',[]);
%set(gca,'ytick',[]);


% Function for painting static data
function stability_paint(handles)
axes(handles.stability);
% f(x)=1/x, splitted to 2 subfunctions
x=str2double(get(handles.editG1,'String'));
y=str2double(get(handles.editG2,'String'));
% stability function
x1=0.333:0.01:3.04;
y1=1./x1;
x2=-3.0:0.01:-0.33;
y2=1./x2;
% filled area
p3=fill([0,0,x1,3],[0,3,y1,0],[0.61 0.61 0.99], ...
        [-3,x2,0,0],[0,y2,-3,0],[0.61 0.61 0.99]);   % background color
set(p3,'EdgeColor','None'); % no border
hold on;
p1=plot(x1,y1,'k', ...              % border
    x2,y2,'k', ...                  % border
    -1:0.01:1,-1:0.01:1,':r', ...   % dotted line
    -3:0.01:3,0,'k', ...            % x-axis
    0,-3:0.01:3,'k', ...            % y-axis
    1,1,'--ro', ...                 
    0,0,'--ro', ...
    0,1,'--ro', ...
    -1,-1,'--ro', ...
    -1,0,'--ro', ...
    2,1/3,'--ro', ...
    'LineWidth',1, ...
    'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
p2=plot(x,y,'--xb','LineWidth',2,'MarkerSize',8);   % pointing cross
hold off;
% labels
text(2.55,-0.25,'g1','Color','k');
text(-0.45,2.75,'g2','Color','k');
% axis constraints
axis([-3 3 -3 3]);
% mouse click handlers (for all graphics elements)
set(gca,'ButtonDownFcn',{@setGs,handles});
set(p1,'ButtonDownFcn',{@setGs,handles});
set(p2,'ButtonDownFcn',{@setGs,handles});
set(p3,'ButtonDownFcn',{@setGs,handles});


% Function to set g1,g2 values by clicking a stability graph
function setGs(hObject, eventdata, handles)
ps=get(gca,'CurrentPoint');
n=1000; % round constant
g1=ps(1);
g2=ps(4);
L=get(handles.sliderL,'value');
set(handles.editG1,'String',num2str(round(n*g1)/n));
set(handles.editG2,'String',num2str(round(n*g2)/n));
set(handles.editR1,'String',num2str(round(n*(L/(1-g1)))/n));
set(handles.editR2,'String',num2str(round(n*(L/(1-g2)))/n));
resonator_paint(handles);
stability_paint(handles);


% Function for painting static data
function energy_paint(handles)
axes(handles.energy);
axis([-1 1 -1 1]);


% --- Outputs from this function are returned to the command line.
function varargout = resonator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function help_menu_Callback(hObject, eventdata, handles)
% hObject    handle to help_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function help_content_Callback(hObject, eventdata, handles)
% hObject    handle to help_content (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Obsah nápovědy','Nápověda');

% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'Fyzikální optika 2012/2013 - Laserový rezonátor'
        ''
        'Autoři: Radek Fér <xferra00@stud.fit.vutbr.cz>'
        '           Miroslav Skácel <xskace00@stud.fit.vutbr.cz>'}, ...
        'O programu');

% --- Executes on sliderL movement.
function sliderL_Callback(hObject, eventdata, handles)
% hObject    handle to sliderL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of sliderL
%        get(hObject,'Min') and get(hObject,'Max') to determine range of sliderL
n = 1000; % round constant
L = get(hObject,'value');
set(handles.editL,'String',num2str(round(n*L)/n));
r1=str2double(get(handles.editR1,'String'));
r2=str2double(get(handles.editR2,'String'));
set(handles.editG1,'String',num2str(round(n*(1-(L/r1)))/n));
set(handles.editG2,'String',num2str(round(n*(1-(L/r2)))/n));
resonator_paint(handles);
stability_paint(handles);

% --- Executes during object creation, after setting all properties.
function sliderL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: sliderL controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function editG1_Callback(hObject, eventdata, handles)
% hObject    handle to editG1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editG1 as text
%        str2double(get(hObject,'String')) returns contents of editG1 as a double
n=1000;
L=get(handles.sliderL,'value');
g1=str2num(get(handles.editG1,'String'));
set(handles.editR1,'String',num2str(round(n*L/(1-g1))/n));
resonator_paint(handles);
stability_paint(handles);
            
% --- Executes during object creation, after setting all properties.
function editG1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editG1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editG2_Callback(hObject, eventdata, handles)
% hObject    handle to editG2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editG2 as text
%        str2double(get(hObject,'String')) returns contents of editG2 as a double
n=1000;
L=get(handles.sliderL,'value');
g2=str2num(get(handles.editG2,'String'));
set(handles.editR2,'String',num2str(round(n*L/(1-g2))/n));
resonator_paint(handles);
stability_paint(handles);

% --- Executes during object creation, after setting all properties.
function editG2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editG2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editL_Callback(hObject, eventdata, handles)
% hObject    handle to editL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editL as text
%        str2double(get(hObject,'String')) returns contents of editL as a double


% --- Executes during object creation, after setting all properties.
function editL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editR1_Callback(hObject, eventdata, handles)
% hObject    handle to editR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editR1 as text
%        str2double(get(hObject,'String')) returns contents of editR1 as a double
n=1000;
L=get(handles.sliderL,'value');
r1=str2num(get(handles.editR1,'String'));
set(handles.editG1,'String',num2str(1-(L/r1)));
resonator_paint(handles);
stability_paint(handles);

% --- Executes during object creation, after setting all properties.
function editR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editR2_Callback(hObject, eventdata, handles)
% hObject    handle to editR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editR2 as text
%        str2double(get(hObject,'String')) returns contents of editR2 as a double
n=1000;
L=get(handles.sliderL,'value');
r2=str2num(get(handles.editR2,'String'));
set(handles.editG2,'String',num2str(1-(L/r2)));
resonator_paint(handles);
stability_paint(handles);

% --- Executes during object creation, after setting all properties.
function editR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


