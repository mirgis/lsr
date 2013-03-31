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

% Last Modified by GUIDE v2.5 30-Mar-2013 23:52:01

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
resonator_init(handles);
stability_init(handles);
energy_init(handles);

% Choose default command line output for resonator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes resonator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% Function for resonator field initialization
function resonator_init(handles)
axes(handles.resonator);
plot(-5:0.1:5,0,':k');
set(gca,'xtick',[]);
set(gca,'ytick',[]);


% Function for stability field initialization
function stability_init(handles)
stability_paint(handles,0.5,0.5);


% Function for painting static data
function stability_paint(handles,x,y)
axes(handles.stability);
% f(x)=1/x, splitted to 2 subfunctions
x1=0.333:0.01:3.0;
x2=-3.0:0.01:-0.33;
y1=1./x1;
y2=1./x2;

 p3=fill([0,0,x1,3],[0,3,y1,0],[0.61 0.61 0.99], ...
         [-3,x2,0,0],[0,y2,-3,0],[0.61 0.61 0.99]);   % background

set(p3,'EdgeColor','None');
hold on;
p1=plot(x1,y1,'k', ...              % points and border lines
    x2,y2,'k', ...
    -1:0.01:1,-1:0.01:1,':r', ...
    -3:0.01:3,0,'k', ...
    0,-3:0.01:3,'k', ...
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

text(2.55,-0.25,'g1','Color','k');
text(-0.45,2.75,'g2','Color','k');
axis([-3 3 -3 3]);
set(gca,'ButtonDownFcn',{@getGs,handles});
set(p1,'ButtonDownFcn',{@getGs,handles});
set(p2,'ButtonDownFcn',{@getGs,handles});
set(p3,'ButtonDownFcn',{@getGs,handles});
set(handles.textL,'String',['L = ' num2str(get(handles.sliderL,'value'))]);


% Function for energy field initialization
function energy_init(handles)
axes(handles.energy);
axis([-1 1 -1 1]);


% Function to get g1,g2 values by clicking a stability graph
function getGs(src,evt,h)
ps=get(gca,'CurrentPoint');
stability_paint(h,ps(1),ps(4));
n=1000;
set(h.textG1,'String',['g1 = ' num2str(round(n*ps(1))/n)]);
set(h.textG2,'String',['g2 = ' num2str(round(n*ps(4))/n)]);


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
n=1000;
set(handles.textL,'String',['L = ' num2str(round(n*get(hObject,'value'))/n)]);

% --- Executes during object creation, after setting all properties.
function sliderL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: sliderL controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on mouse press over axes background.
function stability_Callback(hObject, eventdata, handles)
% hObject    handle to stability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
