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

% Last Modified by GUIDE v2.5 03-Apr-2013 00:24:17

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
global img;
% init
img=imread('stability.png');
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
t1=-PI:0.0001:0; % M1 arc base
t2=0:0.0001:PI;  % M2 arc base
x1=sin(t1);
y1=cos(t1);
x2=sin(t2);
y2=cos(t2);

thr=min(min(abs(r1/3),abs(r2/3)),0.63); % rendering arc threshold

xx=x1((r1*y1>-thr)&(r1*y1<thr)); % M1 x arc trim
yy=y1((r1*y1>-thr)&(r1*y1<thr)); % M1 y arc trim

xxx=x2((r2*y2>-thr)&(r2*y2<thr));% M2 x arc trim 
yyy=y2((r2*y2>-thr)&(r2*y2<thr));% M2 y arc trim

plot(r1*xx+r1-L/2,r1*yy,'b','LineWidth',2,'Color','g'); % M1 mirror
hold on; axis equal;
plot(r1/2-L/2,0,'xg'); % M1 focus
if(r1>0) % M1 concave
 line([r1*xx(1)+r1-L/2, r1-L/2],[r1*(yy(1)), 0],'LineStyle',':');
 line([r1*xx(end)+r1-L/2, r1-L/2],[r1*(yy(end)), 0],'LineStyle',':');
else     % M1 convex
 line([r1*xx(1)+r1-L/2, -r1-L/2],[r1*(yy(1)), 2*r1*(yy(1))],'LineStyle',':');
 line([r1*xx(end)+r1-L/2, -r1-L/2],[r1*(yy(end)), 2*r1*(yy(end))],'LineStyle',':');    
end

plot(r2*xxx-r2+L/2,r2*yyy,'b','LineWidth',2); % M2 mirror
plot(-r2/2+L/2,0,'xb'); % M2 focus
if(r2>0) % M2 concave 
 line([r2*xxx(1)-r2+L/2, -r2+L/2],[r2*(yyy(1)), 0],'LineStyle',':');
 line([r2*xxx(end)-r2+L/2, -r2+L/2],[r2*(yyy(end)), 0],'LineStyle',':');
else     % M2 convex
 line([r2*xxx(1)-r2+L/2, r2+L/2],[r2*(yyy(1)),2*r2*(yyy(1))],'LineStyle',':');
 line([r2*xxx(end)-r2+L/2, r2+L/2],[r2*(yyy(end)), 2*r2*(yyy(end))],'LineStyle',':');   
end
line([-100,100],[0,0],'LineStyle',':','Color','k'); % optical axis
axis([min(-2,min(1.1*(r1*xx+r1-L/2))) max(2,max(1.1*(r2*xxx-r2+L/2))) ... % axis constraint
      min(-0.7,-1.1*thr) max(0.7,1.1*thr)]);
hold off;

% Update stable/unstable
g1 = str2double(get(handles.editG1, 'String'));
g2 = str2double(get(handles.editG2, 'String'));
stable = g1*g2 >= 0 && g1*g2 <= 1;
if stable
    set(handles.text15, 'String', 'Stable');
    set(handles.text15, 'BackgroundColor', [0 1. 0]);
else
    set(handles.text15, 'String', 'Unstable');
    set(handles.text15, 'BackgroundColor', [1. 0 0]);
end

% Graphical stability interpretation using circles with half radii
global waist;
if get(handles.checkbox1, 'Value')
    hold on;
    circle(r1/2-L/2,0,r1/2, 'g');
    circle(-r2/2+L/2,0,r2/2, 'b');
    hold off;

    if stable
        hold on;
        % Compute the intersection of these two circles and draw line
        % at w_0 if there is real solution
        syms X Y R1 R2 X1 X2;
        R1 = r1/2;
        R2 = r2/2;
        X1 = r1/2-L/2;
        Y1 = 0;
        Y2 = 0;
        X2 = -r2/2+L/2;
        circle1 = (X-X1)^2 + (Y-Y1)^2 - R1^2;
        circle2 = (X-X2)^2 + (Y-Y2)^2 - R2^2;
        S=solve([circle1, circle2], 'Real', true);
        if (~isempty(S))
            if length(S.X) == 2
                waist=eval(S.X(1));
                line(eval(S.X), eval(S.Y), 'LineWidth', 2, 'LineStyle', '--');
            end
        end
        hold off;
    end
end

% plot position of z for energy computation
hold on;
z=get(handles.slider3, 'value');
line([z z], [-0.2 0.2]);
hold off;

% Function for painting static data
function stability_paint(handles)
global img;
global waist;

axes(handles.stability);
x=str2double(get(handles.editG1,'String'));
y=str2double(get(handles.editG2,'String'));
i=imagesc([-3 3], [3 -3], img);
hold on; axis xy;
p=plot(x,y,'--xb','LineWidth',2,'MarkerSize',8);   % pointing cross
hold off;
% axis constraints
axis([-3 3 -3 3]);
% mouse click handlers (for all graphics elements)
set(gca,'ButtonDownFcn',{@setGs,handles});
set(p,'ButtonDownFcn',{@setGs,handles});
set(i,'ButtonDownFcn',{@setGs,handles});


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
energy_paint(handles);


% Function for painting static data
function energy_paint(handles)
axes(handles.energy);
global waist;
if isempty(waist)
    waist=0;
end
r=[-0.001:0.000001:0.001];
z=get(handles.slider3, 'value')-waist;
if z == 0
    z = 0.00001;
end
E0=1;
g1 = str2double(get(handles.editG1, 'String'));
g2 = str2double(get(handles.editG2, 'String'));
if ~(0 <= g1*g2 && 1 >= g1*g2)
    cla
    return;
end
L = get(handles.sliderL,'value');
lambda=550*10^-9;
k=2*pi/lambda;
raleyigh_range_sq = L^2*g1*g2*(1-g1*g2)/(g1+g2-2*g1*g2)^2;
R_z = z*(1+(sqrt(raleyigh_range_sq)/z)^2);
w0_sq = (L*lambda/pi)*sqrt(g1*g2*(1-g1*g2)/(g1+g2-2*g1*g2)^2);
w_z=sqrt(w0_sq)*sqrt(1+(z/sqrt(raleyigh_range_sq))^2);
E=E0*(sqrt(w0_sq)/w_z)*exp(-r.^2/w_z^2 - j*k*z - j*k*(r.^2/(2*R_z)) + j*atan(z/sqrt(raleyigh_range_sq)));
plot(r, abs(E));
xlabel('Radial distance [m]')
axis([-0.001 0.001 0 1.2]);
grid;


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
energy_paint(handles);

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
g1=str2double(get(handles.editG1,'String'));
if(g1==1) %
g1=0.998;
set(handles.editG1,'String',g1);
end
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
if(g2==1) %
g2=0.998;
set(handles.editG2,'String',g2);
end
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
n = 1000; % round constant
L = str2double(get(hObject,'String'));
set(handles.sliderL,'value',L);
r1=str2double(get(handles.editR1,'String'));
r2=str2double(get(handles.editR2,'String'));
set(handles.editG1,'String',num2str(round(n*(1-(L/r1)))/n));
set(handles.editG2,'String',num2str(round(n*(1-(L/r2)))/n));
resonator_paint(handles);
stability_paint(handles);

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

function circle(x,y,r, color)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp, 'Color', color, 'LineStyle', '--');


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
resonator_paint(handles);


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
energy_paint(handles);
resonator_paint(handles);


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in anim.
function anim_Callback(hObject, eventdata, handles)
% hObject    handle to anim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop;
axes(handles.resonator);
stop=false;
PI=3.1415926535;

L=get(handles.sliderL,'value'); % mirrors distance
r1=str2double(get(handles.editR1,'String')); % M1 radius
r2=str2double(get(handles.editR2,'String')); % M2 radius
t1=-PI:0.0001:0; % M1 arc base
t2=0:0.0001:PI;  % M2 arc base
x1=sin(t1);
y1=cos(t1);
x2=sin(t2);
y2=cos(t2);

thr=min(min(abs(r1/3),abs(r2/3)),0.63); % rendering arc threshold

xx=x1((r1*y1>-thr)&(r1*y1<thr)); % M1 x arc trim
yy=y1((r1*y1>-thr)&(r1*y1<thr)); % M1 y arc trim

xxx=x2((r2*y2>-thr)&(r2*y2<thr));% M2 x arc trim 
yyy=y2((r2*y2>-thr)&(r2*y2<thr));% M2 y arc trim

x=[0 0 0 0 0 0 0 0];    % start positions
y=[0 0 0 0 0 0 0 0];

dl1=round(length(xx)/3);
dl2=round(length(xxx)/3);

nn1=abs(round(1250/r1));
nn2=abs(round(1250/r2));

xo=[r1*xx(nn1)+r1-L/2 ...
    r1*xx(dl1)+r1-L/2 ...
    r1*xx(2*dl1)+r1-L/2 ...
    r1*xx(end-nn1)+r1-L/2 ...
    r2*xxx(nn2)-r2+L/2 ...
    r2*xxx(dl2)-r2+L/2 ...
    r2*xxx(2*dl2)-r2+L/2 ...    
    r2*xxx(end-nn2)-r2+L/2
    ];   % offset
yo=[r1*(yy(nn1)) ...
    r1*(yy(dl1)) ...
    r1*(yy(2*dl1)) ...
    r1*(yy(end-nn1))...
    r2*(yyy(nn2))...
    r2*(yyy(dl2)) ...
    r2*(yyy(2*dl2))...    
    r2*(yyy(end-nn2))
    ];

phi=[atan((-r1*(yy(nn1)))/(r1-L/2-(r1*xx(nn1)+r1-L/2))) ...
     atan((-r1*(yy(dl1)))/(r1-L/2-(r1*xx(dl1)+r1-L/2)))...
     atan((-r1*(yy(2*dl1)))/(r1-L/2-(r1*xx(dl1)+r1-L/2)))...
     atan((-r1*(yy(end-nn1)))/(r1-L/2-(r1*xx(end-nn1)+r1-L/2))) ...
     PI+atan((-r2*(yyy(nn2)))/(-r2+L/2-(r2*xxx(nn2)-r2+L/2))) ...
     PI+atan((-r2*(yyy(dl2)))/(-r2+L/2-(r2*xxx(dl2)-r2+L/2))) ...
     PI+atan((-r2*(yyy(2*dl2)))/(-r2+L/2-(r2*xxx(2*dl2)-r2+L/2))) ...
     PI+atan((-r2*(yyy(end-nn2)))/(-r2+L/2-(r2*xxx(end-nn2)-r2+L/2)))     
     ];

dx=cos(phi).*.7;
dy=sin(phi).*.7;

step=0.15;
k = [0 0 0 0 0 0 0 0];
refl=[0 0 0 0 0 0 0 0];


while(true)
%for s=1:1:300    
    if(stop)
        break;
    end
plot(r1*xx+r1-L/2,r1*yy,'b','LineWidth',2,'Color','g'); % M1 mirror
hold on; axis equal;
plot(r1/2-L/2,0,'xg'); % M1 focus
if(r1>0) % M1 concave
 line([r1*xx(1)+r1-L/2, r1-L/2],[r1*(yy(1)), 0],'LineStyle',':');
 line([r1*xx(end)+r1-L/2, r1-L/2],[r1*(yy(end)), 0],'LineStyle',':');
else     % M1 convex
 line([r1*xx(1)+r1-L/2, -r1-L/2],[r1*(yy(1)), 2*r1*(yy(1))],'LineStyle',':');
 line([r1*xx(end)+r1-L/2, -r1-L/2],[r1*(yy(end)), 2*r1*(yy(end))],'LineStyle',':');    
end

plot(r2*xxx-r2+L/2,r2*yyy,'b','LineWidth',2); % M2 mirror
plot(-r2/2+L/2,0,'xb'); % M2 focus
if(r2>0) % M2 concave 
 line([r2*xxx(1)-r2+L/2, -r2+L/2],[r2*(yyy(1)), 0],'LineStyle',':');
 line([r2*xxx(end)-r2+L/2, -r2+L/2],[r2*(yyy(end)), 0],'LineStyle',':');
else     % M2 convex
 line([r2*xxx(1)-r2+L/2, r2+L/2],[r2*(yyy(1)),2*r2*(yyy(1))],'LineStyle',':');
 line([r2*xxx(end)-r2+L/2, r2+L/2],[r2*(yyy(end)), 2*r2*(yyy(end))],'LineStyle',':');   
end
line([-100,100],[0,0],'LineStyle',':','Color','k'); % optical axis
  
	dx=cos(phi).*.7;
	dy=sin(phi).*.7;

	x=x+step*dx;
	y=y+step*dy;

    for i=1:length(x)
        if(dx(i)>0)
            Mx = xxx(round(nn2*r2*yyy)/nn2==round(nn2*(yo(i)+y(i)+dy(i)))/nn2);
            refl(i)=inf;
            if(~isempty(Mx))
                refl(i)=L*Mx(1);
            end
            k(i)=-((-(yo(i)+y(i)+dy(i)))/(xo(i)+x(i)+dx(i)+(r2-L/2)));
            quiver(xo(i)+x(i),yo(i)+y(i),dx(i),dy(i));
            if(xo(i)+x(i)+dx(i)+L/2 > refl(i))
                axs=atan(k(i));
                phi(i)=(axs+(PI-phi(i)))+axs;
                x(i)=x(i)+dx(i);
                y(i)=y(i)+dy(i);
            end
        else
            Mx = xx(round(nn1*r1*yy)/nn1==round(nn1*(yo(i)+y(i)+dy(i)))/nn1);
            refl(i)=-inf;
            if(~isempty(Mx))
                refl(i)=(r1-L/2)*Mx(1)+r1-L/2;
            end
            k(i)=-((-(yo(i)+y(i)+dy(i)))/(xo(i)+x(i)+dx(i)-(r1-L/2)));
            quiver(xo(i)+x(i),yo(i)+y(i),dx(i),dy(i));
            if(xo(i)+x(i)+dx(i)+L/2 < refl(i))
                axs=atan(k(i));
                phi(i)=(axs+(PI-phi(i)))+axs;
                x(i)=x(i)+dx(i);
                y(i)=y(i)+dy(i);
            end
        end
    end
   
    axis([min(-2,min(1.1*(r1*xx+r1-L/2))) max(2,max(1.1*(r2*xxx-r2+L/2))) ... % axis constraint
      min(-0.7,-1.1*thr) max(0.7,1.1*thr)]);
    hold off;
    pause(0.01);
end

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop;
stop=true;
