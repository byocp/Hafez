function varargout = proj_ICORD(varargin)
% PROJ_ICORD M-file for proj_ICORD.fig
%      PROJ_ICORD, by itself, creates a new PROJ_ICORD or raises the
%      existing
%      singleton*.
%
%      H = PROJ_ICORD returns the handle to a new PROJ_ICORD or the handle to
%      the existing singleton*.
%
%      PROJ_ICORD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJ_ICORD.M with the given input arguments.
%
%      PROJ_ICORD('Property','Value',...) creates a new PROJ_ICORD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before proj_ICORD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to proj_ICORD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help proj_ICORD

% Last Modified by GUIDE v2.5 20-Oct-2011 12:40:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @proj_ICORD_OpeningFcn, ...
    'gui_OutputFcn',  @proj_ICORD_OutputFcn, ...
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

% --- Executes just before proj_ICORD is made visible.
function proj_ICORD_OpeningFcn(hObject, eventdata, handles, varargin)
global param;
global Flag;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to proj_ICORD (see VARARGIN)

% Choose default command line output for proj_ICORD
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
cla(handles.axes4,'reset');
cla(handles.axes5,'reset');
cla(handles.axes6,'reset');
cla(handles.axes7,'reset');
cla(handles.axes9,'reset');
cla(handles.axes11,'reset');
cnames = {'Mean','Std'};
rnames = {'CVT','CST','AST','BSP'};
set(handles.uitable5,'ColumnName',cnames);
set(handles.uitable5,'RowName',rnames);

set(handles.uitable5,'Data',[]);

Flag.stopcrt = 0;

param.title= [{'C_{a}'},{'R^{0}_{a}'},{'\DeltaV'},{'\tau'},{'V_H'},{'\beta_{H}'},{'\alpha'},{'\gamma'},{'\delta_{H}'},{'H_0'},{'P_{sp}'},{'\alpha_{sp}'}];

axes(handles.axes2)
xlabel('Time [sec]','fontsize',10);
ylabel('Heart Rate [bpm]','fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes4)
xlabel('Time [sec]','fontsize',10);
ylabel('Blood Pressure [mmHg]','fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes3)
xlabel('Time [sec]','fontsize',10);
ylabel(param.title(5),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes5)
xlabel('Time [sec]','fontsize',10);
ylabel(param.title(6),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes6)
xlabel('Time [sec]','fontsize',10);
ylabel(param.title(7),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes7)
xlabel('Time [sec]','fontsize',10);
ylabel(param.title(11),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes9)
xlabel('Time [sec]','fontsize',10);
ylabel('Stability Index','fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes11)
xlabel('Time [sec]','fontsize',10);
grid on
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using proj_ICORD.
if strcmp(get(hObject,'Visible'),'off')
    
end

% UIWAIT makes proj_ICORD wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = proj_ICORD_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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
set(hObject, 'String', {'AC02S','AA02S','AB02S','BH02S','CC02S','CG02S','DW02S','GB02S','GS02S','JF02S','JH02S','MSE02S','MSG02S','RM02S','RS02S','SB02S','SP02S'});


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)

% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)

% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4
% --- Executes on button press in Estimate.
function Estimate_Callback(hObject, eventdata, handles)
%% Estimation
% hObject    handle to Estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% *****************************************************************
% Clearing
% *****************************************************************
cla(handles.axes2,'reset');
cla(handles.axes4,'reset');
cla(handles.axes3,'reset');
cla(handles.axes5,'reset');
cla(handles.axes6,'reset');
cla(handles.axes7,'reset');
cla(handles.axes9,'reset');
cla(handles.axes11,'reset');
set(handles.uitable5,'Data',[]);

% *****************************************************************
% Importing data
% *****************************************************************
global Flag;
data= handles.data;
HR  = cell2mat(data(1));
BP  = cell2mat(data(2));
Fs  = cell2mat(data(3));
LoW = cell2mat(data(4));
dl  = cell2mat(data(5));

t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

BP1=[]; HR1=[]; CO1=[];
ww=LoW*Fs; ov=0;
for i=floor(abs((t1*Fs-ww)/(ww-ov)))+1:floor((t2*Fs-ww)/(ww-ov))+1
    BP1 =[BP1; ones(LoW,1).*mean(BP((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
    HR1 =[HR1; ones(LoW,1).*mean(HR((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
    CO1 =[CO1; ones(LoW,1).*5];
end

handles.data_stairs=[{HR1},{BP1}];
% *****************************************************************
% Initializing
% *****************************************************************
global param;
param.LoW=LoW;
param.Fs=Fs;
param.BP_Msrd_init=BP1(1);
param.HR_Msrd_init=HR1(1);
Flag.stopcrt=0;


bound.ul=3;   %Upper Bound
bound.ll=0.3; %Lower Bound
F1_Initialization({bound.ul,bound.ll});

% *****************************************************************
% Plotting the measured data
% *****************************************************************
axes(handles.axes2)
plot(t1:t1+tend-1,HR1(1:tend),'r','LineWidth',2)
xlim([t1 t2]);
ylim([mean(HR1(1:tend))-3*std(HR1(1:tend)),mean(HR1(1:tend))+3*std(HR1(1:tend))]);
xlabel('Time [sec]','fontsize',10);
ylabel('Heart Rate [bpm]','fontsize',10);
grid on
guidata(hObject, handles); %updates the handles

axes(handles.axes4)
plot(t1:t1+tend-1,BP1(1:tend),'r','LineWidth',2)
xlim([t1 t2]);
ylim([mean(BP1(1:tend))-3*std(BP1(1:tend)),mean(BP1(1:tend))+3*std(BP1(1:tend))]);
xlabel('Time [sec]','fontsize',10);
ylabel('Blood Pressure [mmHg]','fontsize',10);
grid on
guidata(hObject, handles);

% *****************************************************************
% Estimation of the First Batch & other parts
% *****************************************************************
clear HR_W; clear BP_W; clear SV_W; clear N_W; clear pp;
HR1=HR1./60;
SV1 =CO1*1000./(HR1.*60);
[HR_W N_W]=F2_Window(HR1,LoW,0);
BP_W      =F2_Window(BP1,LoW,0);
SV_W      =F2_Window(SV1,LoW,0);
Stab_ind=zeros(1,N_W*LoW);
p(1,1:12)=zeros(1,12);
j=1;
while j<3
    i=1;
    p(i,3)=mean(SV_W(:,i));
    param.ub(3)=p(i,3);    param.lb(3)=p(i,3);
    if j==1
        out=F1_Runopt(BP_W(:,i),HR_W(:,i),param.nom);
    end
    
    if j==2
        out=F1_Runopt(BP_W(:,i),HR_W(:,i),p(i,:));
    end
    p(i,:)=cell2mat(out(1));
    Q=F1_Sim([p(1,:),param.BP_Msrd_init,param.HR_Msrd_init]);
    x=cell2mat(Q(1));
    y=cell2mat(Q(2));
    param.BP_Msrd_init=x(end);
    param.HR_Msrd_init=y(end);
    x_sim(LoW*(i-1)+1:LoW*i,1) = x(30-LoW+1:30);
    y_sim(LoW*(i-1)+1:LoW*i,1) = y(30-LoW+1:30);
    j=j+1;
end

pp(5,1:LoW) = p(1,5)      * ones(1,LoW).*1./(1+exp(-p(1,12).*(BP1(1+LoW*(1-1))-p(1,11))));
pp(6,1:LoW) = p(1,6)      * ones(1,LoW).*(1-1./(1+exp(-p(1,12).*(BP1(1+LoW*(1-1))-p(1,11)))));
pp(7,1:LoW) = p(1,7)      * ones(1,LoW).*(1-1./(1+exp(-p(1,12).*(BP1(1+LoW*(1-1))-p(1,11)))));
pp(11,1:LoW)= p(1,11)     * ones(1,LoW);
Stab_ind    = SI({p(1,:)})* ones(1,LoW);

% *****************************************************************
% Plotting the first batch of simulated data (HR and BP)
% *****************************************************************
axes(handles.axes2)
hold on
plot(t1:t1+LoW-1,y_sim(1:LoW)*60,'k','LineWidth',1);
box on;
guidata(hObject, handles); %updates the handles

axes(handles.axes4)
hold on
plot(t1:t1+LoW-1,x_sim(1:LoW),'k','LineWidth',1);
box on;
guidata(hObject, handles);

for i=2:N_W
    if Flag.stopcrt~=1
        p(i,3)=mean(SV_W(:,i));
        %           param.ub(3)=p(i,3)*bound.ul;   param.lb(3)=p(i,3)*bound.ll;
        param.ub(3)=p(i,3);    param.lb(3)=p(i,3);
        out=F1_Runopt(BP_W(:,i),HR_W(:,i),p(i-1,:));
        p(i,:)=cell2mat(out(1));
        
        param.BP_Msrd_init=x(end);
        param.HR_Msrd_init=y(end);
        
        Q=F1_Sim([p(i,:),param.BP_Msrd_init,param.HR_Msrd_init]);
        
        x=cell2mat(Q(1));
        y=cell2mat(Q(2));
        clear Q;
        x_sim(LoW*(i-1)+1:LoW*i,1) = x(30-LoW+1:30);
        y_sim(LoW*(i-1)+1:LoW*i,1) = y(30-LoW+1:30);
        
        axes(handles.axes2)
        hold on
        plot(t1:t1+LoW*i-1,y_sim(1:LoW*i)*60,'k','LineWidth',1);
        guidata(hObject, handles); %updates the handles
        
        axes(handles.axes4)
        hold on
        plot(t1:t1+LoW*i-1,x_sim(1:LoW*i),'k','LineWidth',1);
        guidata(hObject, handles); %updates the handles
        %**********************************************************************
        pp(5,LoW*(i-1):LoW*i-1)      = p(i,5)*ones(1,LoW).*1./(1+exp(-p(i,12).*(BP1(1+LoW*(i-1))-p(i,11))));
        pp(6,LoW*(i-1):LoW*i-1)      = p(i,6)*ones(1,LoW).*(1-1./(1+exp(-p(i,12).*(BP1(1+LoW*(i-1))-p(i,11)))));
        pp(7,LoW*(i-1):LoW*i-1)      = p(i,7)*ones(1,LoW).*(1-1./(1+exp(-p(i,12).*(BP1(1+LoW*(i-1))-p(i,11)))));
        pp(11,LoW*(i-1):LoW*i-1)     = p(i,11)*ones(1,LoW);
        Stab_ind(1,LoW*(i-1):LoW*i-1)= SI({p(i,:)})*ones(1,LoW);
        %******************************************************************
        axes(handles.axes3)
        xlim([t1 t2]);
        hold on
        plot(t1:t1+size(pp,2)-1,pp(5,:),'k');
        grid on;
        box on;
        xlabel('Time [sec]','fontsize',10);
        ylabel(param.title(5),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes5)
        xlim([t1 t2]);
        hold on
        plot(t1:t1+size(pp,2)-1,pp(6,:),'k');
        grid on;
        box on;
        xlabel('Time [sec]','fontsize',10);
        ylabel(param.title(6),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes6)
        xlim([t1 t2]);
        hold on;
        plot(t1:t1+size(pp,2)-1,pp(7,:),'k');
        grid on;
        box on;
        xlabel('Time [sec]','fontsize',10);
        ylabel(param.title(7),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes7)
        xlim([t1 t2]);
        hold on;
        plot(t1:t1+size(pp,2)-1,pp(11,:),'k');
        grid on;
        box on;
        xlabel('Time [sec]','fontsize',10);
        ylabel(param.title(11),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes11)
        xlim([t1 t2]);
        plot(t1:t1+size(pp,2)-1,pp(7,:),'r');
        hold on
        plot(t1:t1+size(pp,2)-1,pp(6,:),'-b');
        hold on
        plot(t1:t1+size(pp,2)-1,pp(5,:),'.k');
        hold on
        grid on;
        box on;
        xlabel('Time [sec]','fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes9)
        xlim([t1 t2]);
        hold on;
        plot(t1:t1+size(Stab_ind,2)-1,Stab_ind,'k');
        grid on;
        box on
        xlabel('Time [sec]','fontsize',10);
        ylabel('Stability Index','fontsize',10);
        guidata(hObject, handles); %updates the handles
        %**********************************************************************
        NumDat=[mean(pp(5,:)),std(pp(5,:)); mean(pp(6,:)),std(pp(6,:)); ...
            mean(pp(7,:)),std(pp(7,:)); mean(pp(11,:)),std(pp(11,:));];
        set(handles.uitable5,'Data',NumDat);
        uitable5_CellEditCallback(hObject, eventdata, handles)
    end
end

%**********************************************************************
axes(handles.axes2)
HRres_legend=legend('Measurement','Simulation');
set(HRres_legend,'FontSize',8);
guidata(hObject, handles); %updates the handles
%**********************************************************************
axes(handles.axes4)
BPres_legend=legend('Measurement','Simulation');
set(BPres_legend,'FontSize',8);
guidata(hObject, handles);

handles.Simulation=[{x_sim},{y_sim}];
handles.Estimation=pp;

guidata(hObject, handles);


%**********************************************************************
%Just in case of making smooth data
% axes(handles.axes3)
% pp5_smooth=F2_Smooth([{pp(5,:)'},{2},{0},{0}]);
% hold on
% plot(pp5_smooth,'r');
% guidata(hObject, handles); %updates the handles

% --- Executes on button press in Plot.
function Plot_Callback(hObject, eventdata, handles)
cla(handles.axes2,'reset');
cla(handles.axes4,'reset');

% Importing data
data=handles.data;
HR=cell2mat(data(1));
BP=cell2mat(data(2));
Fs=cell2mat(data(3));
LoW=cell2mat(data(4));
dl=cell2mat(data(5));

t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

BP1=[]; HR1=[]; CO1=[];

ww=LoW*Fs; ov=0;
for i=floor(abs((t1*Fs-ww)/(ww-ov)))+1:floor((t2*Fs-ww)/(ww-ov))+1
    BP1=[BP1; ones(LoW,1).*mean(BP((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
    HR1=[HR1; ones(LoW,1).*mean(HR((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
    CO1=[CO1; ones(LoW,1).*5];
end


%% Plotting data
axes(handles.axes2)
plot(t1:1/Fs:t2,HR(t1*Fs:t2*Fs));
hold on
plot(t1:t1+LoW*floor((t2-t1)/LoW)-1,HR1(1:LoW*floor((t2-t1)/LoW)),'r','LineWidth',2);
ylim([mean(HR1(1:tend))-3*std(HR1(1:tend)),mean(HR1(1:tend))+3*std(HR1(1:tend))]);
xlim([t1 t2]);
xlabel('Time [sec]','fontsize',10);
ylabel('Heart Rate [bpm]','fontsize',10);
HR_legend=legend('Instantaneous','Mean');
set(HR_legend,'FontSize',8);
grid on
guidata(hObject, handles); %updates the handles

axes(handles.axes4)
plot(t1:1/Fs:t2,BP(t1*Fs:t2*Fs));
hold on
plot(t1:t1+LoW*floor((t2-t1)/LoW)-1,BP1(1:LoW*floor((t2-t1)/LoW)),'r','LineWidth',2);
ylim([mean(BP1(1:tend))-3*std(BP1(1:tend)),mean(BP1(1:tend))+3*std(BP1(1:tend))]);
xlim([t1 t2]);
xlabel('Time [sec]','fontsize',10);
ylabel('Blood Pressure [mmHg]','fontsize',10);
BP_legend=legend('Instantaneous','Mean');
set(BP_legend,'FontSize',8);
grid on
guidata(hObject, handles);
%updates the handles
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sliderValue = get(handles.slider1,'Value');
set(handles.slider1_editText,'String', num2str(sliderValue));
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider1_editText_Callback(hObject, eventdata, handles)
% hObject    handle to slider1_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sliderValue = get(handles.slider1_editText,'String');
sliderValue = str2num(sliderValue);
if (isempty(sliderValue) || sliderValue < 0 )
    set(handles.slider1,'Value',0);
    set(handles.slider1_editText,'String','0');
else
    set(handles.slider1,'Value',sliderValue);
end
% Hints: get(hObject,'String') returns contents of slider1_editText as text
%        str2double(get(hObject,'String')) returns contents of slider1_editText as a double


% --- Executes during object creation, after setting all properties.
function slider1_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=handles.data;
Fs=cell2mat(data(3));
dl=cell2mat(data(5));
set(handles.slider2,'Max', floor(dl/Fs));

sliderValue = get(handles.slider2,'Value');
set(handles.slider2_editText,'String',int2str(sliderValue));
% set(handles.slider2_editText,'String', num2str(sliderValue));
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes during object creation, after setting all properties.
function slider2_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function slider2_editText_Callback(hObject, eventdata, handles)
% hObject    handle to slider2_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sliderValue = get(handles.slider2_editText,'String');
sliderValue = str2num(sliderValue);
if (isempty(sliderValue) || sliderValue < 0 )
    set(handles.slider2,'Value',0);
    set(handles.slider2_editText,'String','0');
else
    set(handles.slider2,'Value',sliderValue);
end
% Hints: get(hObject,'String') returns contents of slider2_editText as text
%        str2double(get(hObject,'String')) returns contents of slider2_editText as a double


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


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
set(hObject, 'String', {'1Hz','100Hz','1000Hz'});
guidata(hObject, handles);


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


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
set(hObject, 'String', {'1','5','10','30'});
guidata(hObject, handles);

% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12
popup_sel_index12 = get(handles.popupmenu12, 'Value');
switch popup_sel_index12
    case 1
        set(handles.popupmenu1,'String',{'AC02S','AA02S','AB02S','BH02S','CC02S','CG02S','DW02S','GB02S','GS02S','JF02S','JH02S','MSE02S','MSG02S','RM02S','RS02S','SB02S','SP02S'});
    case 2
        set(handles.popupmenu1,'String',{'486','484','477','476','474','289'});
    case 3
        set(handles.popupmenu1,'String',{});
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Import.
function Import_Callback(hObject, eventdata, handles)
% hObject    handle to Import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
cla(handles.axes4,'reset');
cla(handles.axes5,'reset');
cla(handles.axes6,'reset');
cla(handles.axes7,'reset');
cla(handles.axes9,'reset');
cla(handles.axes11,'reset');
set(handles.uitable5,'Data',[]);

global param;
axes(handles.axes2)
xlabel('Time [sec]','fontsize',10);
ylabel('Heart Rate [bpm]','fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes4)
xlabel('Time [sec]','fontsize',10);
ylabel('Blood Pressure [mmHg]','fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes3)
xlabel('Time [sec]','fontsize',10);
ylabel(param.title(5),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes5)
xlabel('Time [sec]','fontsize',10);
ylabel(param.title(6),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes6)
xlabel('Time [sec]','fontsize',10);
ylabel(param.title(7),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes7)
xlabel('Time [sec]','fontsize',10);
ylabel(param.title(11),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes9)
xlabel('Time [sec]','fontsize',10);
ylabel('Stability Index','fontsize',10);
grid on
guidata(hObject, handles);


axes(handles.axes11)
xlabel('Time [sec]','fontsize',10);
grid on
guidata(hObject, handles);

popup_sel_index12 = get(handles.popupmenu12, 'Value');
popup_sel_index = get(handles.popupmenu1, 'Value');
if popup_sel_index12==1
    switch popup_sel_index
        case 1
            file_name={'AC02S'};
        case 2
            file_name={'AA02S'};
        case 3
            file_name={'AB02S'};
        case 4
            file_name={'BH02S'};
        case 5
            file_name={'CC02S'};
        case 6
            file_name={'CG02S'};
        case 7
            file_name={'DW02S'};
        case 8
            file_name={'GB02S'};
        case 9
            file_name={'GS02S'};
        case 10
            file_name={'JF02S'};
        case 11
            file_name={'JH02S'};
        case 12
            file_name={'MSE02S'};
        case 13
            file_name={'MSG02S'};
        case 14
            file_name={'RM02S'};
        case 15
            file_name={'RS02S'};
        case 16
            file_name={'SB02S'};
        case 17
            file_name={'SP02S'};
    end
    rawdata=importdata(['../../Dataset/ICORD/data/Rugby_games_HRV/dat02_cleaned/',cell2mat(file_name),'.dat']);
    x=rawdata.data;
    dl=length(x);
    if strcmp(file_name,'MSE02S')==1 || strcmp(file_name,'RM02S')==1 || strcmp(file_name,'CC02S')==1 || strcmp(file_name,'GS02S')==1
        BP     =x(1:dl,4);
        HR     =x(1:dl,3);
    else
        BP     =x(1:dl,3);
        HR     =x(1:dl,4);
    end
end

if popup_sel_index12==2
    switch popup_sel_index
        case 1
            file_name={'486'};
        case 2
            file_name={'484'};
        case 3
            file_name={'477'};
        case 4
            file_name={'476'};
        case 5
            file_name={'474'};
        case 6
            file_name={'289'};
    end
    plotATM1=['../../Dataset/MIMIC/',cell2mat(file_name),'nm.mat'];
    plotATM2=['../../Dataset/MIMIC/',cell2mat(file_name),'nm.info'];
    x=plotATM(plotATM1,plotATM2);
    dl=length(x);
    BP     =x(1:dl,1);
    BPsys  =x(1:dl,2);
    BPdias =x(1:dl,3);
    CO     =x(1:dl,5);
    HR     =x(1:dl,8);
    
end

% Remove NaN from signals
HR(isnan(HR)) = 250;
BP(isnan(BP)) = 250;

popup_sel_index5 = get(handles.popupmenu5, 'Value');
switch popup_sel_index5
    case 1
        LoW=1;
    case 2
        LoW=5;
    case 3
        LoW=10;
    case 4
        LoW=30;
end

popup_sel_index3 = get(handles.popupmenu3, 'Value');
switch popup_sel_index3
    case 1
        Fs=1;
    case 2
        Fs=100;
    case 3
        Fs=1000;
end

dl=length(HR);

set(handles.slider1,'Max', round(min(length(HR),length(BP))/Fs));
set(handles.slider1,'Value',0);
set(handles.slider1_editText,'String','0');
set(handles.slider2,'Max', round(min(length(HR),length(BP))/Fs));
set(handles.slider2,'Value',round(min(length(HR),length(BP))/Fs));
set(handles.slider2_editText,'String',int2str(round(min(length(HR),length(BP))/Fs)));

handles.data=[{HR},{BP},{Fs},{LoW},{dl}];
guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in uitable5.
function uitable5_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable5 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



function out=SI(vargin)
C_a       = 1.55;
R0_a      = 0.6;
DeltaV    = 50;
IHR       = 1.66;
tau       = 3;
V_H       = 1.17;
beta_H    = 0.84;
P_init    = 160;
HR_init   = 2;
Alpha     = 1.3;
gamma     = 0.2;
Delta_h   = 1.7;
sig_sp    = 100;
sig_Alpha = 0.05;

p=cell2mat(vargin(1));

input2(1,1)=C_a;
input2(2,1)=R0_a;
input2(3,1)=50;
input2(4,1)=IHR;
input2(5,1)=tau;
input2(6,1)=p(1,5);
input2(7,1)=p(1,6);
input2(8,1)=P_init;
input2(9,1)=HR_init;
input2(10,1)=p(1,7);
input2(11,1)=gamma;
input2(12,1)=Delta_h;
input2(13,1)=p(1,11);
input2(14,1)=sig_Alpha;

out2=F1_Fsolve_Sol(input2);
Pf        = cell2mat(out2(1));
Hf        = cell2mat(out2(2));
outfsolve = cell2mat(out2(3));
out3=F1_Lyap_EigVal([input2;Pf]);
rleigvec=cell2mat(out3(1));
out=max(rleigvec);


% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
global Flag;
% hObject    handle to Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Flag.stopcrt=1;


% --- Executes during object creation, after setting all properties.
function text7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('ICORDPoster_Final_LowRes.pdf')


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
global param;
pp=handles.Estimation;
data= handles.data;
HR  = cell2mat(data(1));
BP  = cell2mat(data(2));
Fs  = cell2mat(data(3));
LoW = cell2mat(data(4));
dl  = cell2mat(data(5));
t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    cla(handles.axes3,'reset');
    axes(handles.axes3)
    hist(pp(5,:),20)
    xlabel(param.title(5),'fontsize',10);
    ylabel('Frequency','fontsize',10);
    guidata(hObject, handles); %updates the handles
    grid on;
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes3,'reset');
    axes(handles.axes3)
    xlim([t1 t2]);
    hold on
    plot(t1:t1+size(pp,2)-1,pp(5,:),'k');
    grid on;
    box on;
    xlabel('Time [sec]','fontsize',10);
    ylabel(param.title(5),'fontsize',10);
    guidata(hObject, handles); %updates the handles
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
global param;
pp=handles.Estimation;
data= handles.data;
HR  = cell2mat(data(1));
BP  = cell2mat(data(2));
Fs  = cell2mat(data(3));
LoW = cell2mat(data(4));
dl  = cell2mat(data(5));
t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    cla(handles.axes5,'reset');
    axes(handles.axes5)
    hist(pp(6,:),20)
    xlabel(param.title(6),'fontsize',10);
    ylabel('Frequency','fontsize',10);
    guidata(hObject, handles); %updates the handles
    grid on;
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes5,'reset');
    axes(handles.axes5)
    xlim([t1 t2]);
    hold on
    plot(t1:t1+size(pp,2)-1,pp(6,:),'k');
    grid on;
    box on;
    xlabel('Time [sec]','fontsize',10);
    ylabel(param.title(6),'fontsize',10);
    guidata(hObject, handles); %updates the handles
    
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
global param;
pp=handles.Estimation;
data= handles.data;
HR  = cell2mat(data(1));
BP  = cell2mat(data(2));
Fs  = cell2mat(data(3));
LoW = cell2mat(data(4));
dl  = cell2mat(data(5));
t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    cla(handles.axes6,'reset');
    axes(handles.axes6)
    hist(pp(7,:),20)
    xlabel(param.title(7),'fontsize',10);
    ylabel('Frequency','fontsize',10);
    guidata(hObject, handles); %updates the handles
    grid on;
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes6,'reset');
    axes(handles.axes6)
    xlim([t1 t2]);
    hold on
    plot(t1:t1+size(pp,2)-1,pp(7,:),'k');
    grid on;
    box on;
    xlabel('Time [sec]','fontsize',10);
    ylabel(param.title(7),'fontsize',10);
    guidata(hObject, handles); %updates the handles
    
end

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
global param;
pp=handles.Estimation;
data= handles.data;
HR  = cell2mat(data(1));
BP  = cell2mat(data(2));
Fs  = cell2mat(data(3));
LoW = cell2mat(data(4));
dl  = cell2mat(data(5));
t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    cla(handles.axes7,'reset');
    axes(handles.axes7)
    hist(pp(11,:),20)
    xlabel(param.title(11),'fontsize',10);
    ylabel('Frequency','fontsize',10);
    guidata(hObject, handles); %updates the handles
    grid on;
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes7,'reset');
    axes(handles.axes7)
    xlim([t1 t2]);
    hold on
    plot(t1:t1+size(pp,2)-1,pp(11,:),'k');
    grid on;
    box on;
    xlabel('Time [sec]','fontsize',10);
    ylabel(param.title(11),'fontsize',10);
    guidata(hObject, handles); %updates the handles
    
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
est=handles.Simulation;
BP_Sim=cell2mat(est(1));
HR_Sim=cell2mat(est(2))*60;

data_stairs= handles.data_stairs;
HR1  = cell2mat(data_stairs(1));
BP1  = cell2mat(data_stairs(2));

data= handles.data;
HR  = cell2mat(data(1));
BP  = cell2mat(data(2));
Fs  = cell2mat(data(3));
LoW = cell2mat(data(4));
dl  = cell2mat(data(5));
t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

err_BP=(BP_Sim-BP1(t1:t1+length(BP_Sim)-1,1))./BP1(t1:t1+length(BP_Sim)-1,1);
err_HR=(HR_Sim-HR1(t1:t1+length(HR_Sim)-1,1))./HR1(t1:t1+length(HR_Sim)-1,1);

if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    cla(handles.axes2,'reset');
    axes(handles.axes2)
    hist(err_HR,20)
    xlabel('Error','fontsize',10);
    ylabel('Frequency','fontsize',10);
    title('Heart Rate');
    grid on;
    box on;
    guidata(hObject, handles); %updates the handles
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes2,'reset');
    axes(handles.axes2)
    plot(t1:t1+tend-1,HR1(1:tend),'r','LineWidth',2)
    hold on
    plot(t1:t1+length(HR_Sim)-1,HR_Sim,'k','LineWidth',1);
    xlim([t1 t2]);
    ylim([mean(HR1(1:tend))-3*std(HR1(1:tend)),mean(HR1(1:tend))+3*std(HR1(1:tend))]);
    xlabel('Time [sec]','fontsize',10);
    ylabel('Heart Rate [bpm]','fontsize',10);
    HRres_legend=legend('Measurement','Simulation');
    set(HRres_legend,'FontSize',8);
    grid on
    box on;
    guidata(hObject, handles); %updates the handles
end

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
est=handles.Simulation;
BP_Sim=cell2mat(est(1));
HR_Sim=cell2mat(est(2));

data_stairs= handles.data_stairs;
HR1  = cell2mat(data_stairs(1));
BP1  = cell2mat(data_stairs(2));

data= handles.data;
HR  = cell2mat(data(1));
BP  = cell2mat(data(2));
Fs  = cell2mat(data(3));
LoW = cell2mat(data(4));
dl  = cell2mat(data(5));
t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

err_BP=(BP_Sim-BP1(t1:t1+length(BP_Sim)-1,1))./BP1(t1:t1+length(BP_Sim)-1,1);
err_HR=(HR_Sim-HR1(t1:t1+length(BP_Sim)-1,1))./HR1(t1:t1+length(BP_Sim)-1,1);

if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    cla(handles.axes4,'reset');
    axes(handles.axes4)
    hist(err_BP,20)
    xlabel('Error','fontsize',10);
    ylabel('Frequency','fontsize',10);
    title('Blood Pressure');
    grid on;
    box on;
    guidata(hObject, handles); %updates the handles
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes4,'reset');
    axes(handles.axes4)
    plot(t1:t1+tend-1,BP1(1:tend),'r','LineWidth',2)
    hold on
    plot(t1:t1+length(BP_Sim)-1,BP_Sim,'k','LineWidth',1);
    xlim([t1 t2]);
    ylim([mean(BP1(1:tend))-3*std(BP1(1:tend)),mean(BP1(1:tend))+3*std(BP1(1:tend))]);
    xlabel('Time [sec]','fontsize',10);
    ylabel('Blood Pressure [mmHg]','fontsize',10);
    BPres_legend=legend('Measurement','Simulation');
    set(BPres_legend,'FontSize',8);
    grid on
    box on;
    guidata(hObject, handles); %updates the handles
end
% Hint: get(hObject,'Value') returns toggle state of checkbox7
