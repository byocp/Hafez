function varargout = ANS_MNT(varargin)
% ANS_MNT M-file for ANS_MNT.fig
%      ANS_MNT, by itself, creates a new ANS_MNT or raises the
%      existing
%      singleton*.
%
%      H = ANS_MNT returns the handle to a new ANS_MNT or the handle to
%      the existing singleton*.
%
%      ANS_MNT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANS_MNT.M with the given input arguments.
%
%      ANS_MNT('Property','Value',...) creates a new ANS_MNT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ANS_MNT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ANS_MNT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ANS_MNT

% Last Modified by GUIDE v2.5 01-Nov-2011 16:09:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ANS_MNT_OpeningFcn, ...
    'gui_OutputFcn',  @ANS_MNT_OutputFcn, ...
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

% --- Executes just before ANS_MNT is made visible.
function ANS_MNT_OpeningFcn(hObject, eventdata, handles, varargin)
global param;
global param_mod;
global Flag;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ANS_MNT (see VARARGIN)

% Choose default command line output for ANS_MNT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
cla(handles.axes4,'reset');
cla(handles.axes5,'reset');
cla(handles.axes6,'reset');
cla(handles.axes11,'reset');
cla(handles.axes12,'reset');
cla(handles.axes13,'reset');
cla(handles.axes14,'reset');
cnames = {'Mean','Std'};
rnames = {'CVT','CST','AST'};
set(handles.uitable5,'ColumnName',cnames);
set(handles.uitable5,'RowName',rnames);

set(handles.uitable5,'Data',[]);

Flag.stopcrt = 0;

param.title= [{'C_{a}'},{'R^{0}_{a}'},{'\DeltaV'},{'\tau'},{'V_H'},{'\beta_{H}'},{'\alpha'},{'\gamma'},{'\delta_{H}'},{'H_0'},{'P_{sp}'},{'\alpha_{sp}'}];
param_mod.title=[{'V_H.T_p'},{'\beta_H.T_s'},'\alpha.T_s'];
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
ylabel(param_mod.title(1),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes5)
xlabel('Time [sec]','fontsize',10);
ylabel(param_mod.title(2),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes6)
xlabel('Time [sec]','fontsize',10);
ylabel(param_mod.title(3),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes11)
xlabel('Time [sec]','fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes12)
xlabel(param_mod.title(1),'fontsize',10);
ylabel(param_mod.title(2),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes13)
xlabel(param_mod.title(1),'fontsize',10);
ylabel(param_mod.title(3),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes14)
xlabel(param_mod.title(2),'fontsize',10);
ylabel(param_mod.title(3),'fontsize',10);
grid on
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using ANS_MNT.
if strcmp(get(hObject,'Visible'),'off')
    
end

% UIWAIT makes ANS_MNT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ANS_MNT_OutputFcn(hObject, eventdata, handles)
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
set(hObject, 'String', {'AC02S','AA02S','AB02S','BH02S','CC02S','CG02S','DW02S','GB02S','GS02S','JF02S','JH02S','MSE02S','MSG02S','RM02S','RS02S','SB02S','SP02S','TH02S','FL02S','TC02S'});


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
cla(handles.axes11,'reset');
cla(handles.axes12,'reset');
cla(handles.axes13,'reset');
cla(handles.axes14,'reset');
set(handles.uitable5,'Data',[]);

% *****************************************************************
% Importing data
% *****************************************************************
global Flag;
global Imp;

data= handles.data;
HR  = cell2mat(data(1));
BP  = cell2mat(data(2));
Fs  = cell2mat(data(3));
LoW = cell2mat(data(4));
dl  = cell2mat(data(5));
if length(data)==6
    CO  = cell2mat(data(6));
end

if length(data)==6
    out=F2_Clean([{HR},{BP},{CO}]);
    HR=cell2mat(out(1));
    BP=cell2mat(out(2));
    CO=cell2mat(out(3));
end

if length(data)==5
    out=F2_Clean([{HR},{BP}]);
    HR=cell2mat(out(1));
    BP=cell2mat(out(2));
end

t1  = max(1,str2double(get(handles.slider1_editText,'String')));
t2  = min(floor(dl/Fs),str2double(get(handles.slider2_editText,'String')));
tend=LoW*floor((t2-t1)/LoW);

BP1=[]; HR1=[]; CO1=[];
ww=LoW*Fs; ov=0;

%without cardiac output measurement
if length(data)==5
    for i=floor(abs((t1*Fs-ww)/(ww-ov)))+1:floor((t2*Fs-ww)/(ww-ov))+1
        BP1 =[BP1; ones(LoW,1).*mean(BP((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
        HR1 =[HR1; ones(LoW,1).*mean(HR((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
        CO1 =[CO1; ones(LoW,1).*5];
    end
end

% with cardiac output measurement
if length(data)==6
    for i=floor(abs((t1*Fs-ww)/(ww-ov)))+1:floor((t2*Fs-ww)/(ww-ov))+1
        BP1 =[BP1; ones(LoW,1).*mean(BP((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
        HR1 =[HR1; ones(LoW,1).*mean(HR((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
        CO1 =[CO1; ones(LoW,1).*mean(CO((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1))];
    end
end

handles.data_stairs=[{HR1},{BP1}];

% *****************************************************************
% Initializing
% *****************************************************************
global param;
global param_mod;
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
param.N_W=N_W;
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
Stab_ind    = F1_SI({p(1,:)})* ones(1,LoW);

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
dist56=zeros(1,N_W);
dist57=zeros(1,N_W);
dist67=zeros(1,N_W);
X56=[];
X57=[];
X67=[];
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
        Stab_ind(1,LoW*(i-1):LoW*i-1)= F1_SI({p(i,:)})*ones(1,LoW);
        %******************************************************************
        axes(handles.axes3)
        xlim([t1 t2]);
        hold on
        plot(t1:t1+size(pp,2)-1,pp(5,:),'c','LineWidth',2);
        grid on;
        box on;
        xlabel('Time [sec]','fontsize',10);
        ylabel(param_mod.title(1),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes5)
        xlim([t1 t2]);
        hold on
        plot(t1:t1+size(pp,2)-1,pp(6,:),'g','LineWidth',2);
        grid on;
        box on;
        xlabel('Time [sec]','fontsize',10);
        ylabel(param_mod.title(2),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes6)
        xlim([t1 t2]);
        hold on;
        plot(t1:t1+size(pp,2)-1,pp(7,:),'b','LineWidth',2);
        grid on;
        box on;
        xlabel('Time [sec]','fontsize',10);
        ylabel(param_mod.title(3),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes11)
        YY(i,:)=[pp(7,LoW*(i-1)-1),pp(6,LoW*(i-1)-1),pp(5,LoW*(i-1)-1)];
        harea=area(YY);
        set(harea(3),'FaceColor','c');
        grid on;
        box on;
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes12)
        scatter(pp(5,LoW*(i-1)-1),pp(6,LoW*(i-1)-1),'k');
        dist56(i)=sqrt((pp(5,LoW*(i-1)-1)-mean(pp(5,:)))^2+(pp(6,LoW*(i-1)-1)-mean(pp(6,:)))^2);
        X56(i-1,:)=[pp(5,LoW*(i-1)-1),pp(6,LoW*(i-1)-1)];
        Z5(i-1)=pp(5,LoW*(i-1)-1);
        hold on
        grid on;
        box on;
        xlabel(param_mod.title(1),'fontsize',10);
        ylabel(param_mod.title(2),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes13)
        scatter(pp(5,LoW*(i-1)-1),pp(7,LoW*(i-1)-1),'k');
        dist57(i)=sqrt((pp(5,LoW*(i-1)-1)-mean(pp(5,:)))^2+(pp(7,LoW*(i-1)-1)-mean(pp(7,:)))^2);
        X57(i-1,:)=[pp(5,LoW*(i-1)-1),pp(7,LoW*(i-1)-1)];
        Z6(i-1)=pp(6,LoW*(i-1)-1);
        hold on
        grid on;
        box on;
        xlabel(param_mod.title(1),'fontsize',10);
        ylabel(param_mod.title(3),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        axes(handles.axes14)
        scatter(pp(6,LoW*(i-1)-1),pp(7,LoW*(i-1)-1),'k');
        dist67(i)=sqrt((pp(6,LoW*(i-1)-1)-mean(pp(6,:)))^2+(pp(7,LoW*(i-1)-1)-mean(pp(7,:)))^2);
        X67(i-1,:)=[pp(6,LoW*(i-1)-1),pp(7,LoW*(i-1)-1)];
        Z7(i-1)=pp(7,LoW*(i-1)-1);
        hold on;
        grid on;
        box on;
        xlabel(param_mod.title(2),'fontsize',10);
        ylabel(param_mod.title(3),'fontsize',10);
        guidata(hObject, handles); %updates the handles
        %******************************************************************
        NumDat=[mean(pp(5,:)),std(pp(5,:)); mean(pp(6,:)),std(pp(6,:)); ...
            mean(pp(7,:)),std(pp(7,:));];
        set(handles.uitable5,'Data',NumDat);
        uitable5_CellEditCallback(hObject, eventdata, handles)
    end
end
Imp.P56{8}=X56;
Imp.P57{8}=X57;
Imp.P67{8}=X67;
Imp.CVT{8}=Z5';
Imp.CST{8}=Z6';
Imp.AST{8}=Z7';

%**************************************************************************
% opts = statset('Display','final');
% [idx,ctrs] = kmeans(X56,2,...
%     'Distance','city',...
%     'Replicates',5,...
%     'Options',opts);
% axes(handles.axes12)
% scatter(X56(idx==1,1),X56(idx==1,2),'r');
% hold on
% scatter(X56(idx==2,1),X56(idx==2,2),'b');
% plot(ctrs(:,1),ctrs(:,2),'kx',...
%     'MarkerSize',8,'LineWidth',2)
% plot(ctrs(:,1),ctrs(:,2),'ko',...
%     'MarkerSize',8,'LineWidth',2)
% xlabel(param_mod.title(1),'fontsize',10);
% ylabel(param_mod.title(2),'fontsize',10);
% guidata(hObject, handles); %updates the handles
% 
% opts = statset('Display','final');
% [idx,ctrs] = kmeans(X57,2,...
%     'Distance','city',...
%     'Replicates',5,...
%     'Options',opts);
% axes(handles.axes13)
% scatter(X57(idx==1,1),X57(idx==1,2),'r');
% hold on
% scatter(X57(idx==2,1),X57(idx==2,2),'b');
% plot(ctrs(:,1),ctrs(:,2),'kx',...
%     'MarkerSize',8,'LineWidth',2)
% plot(ctrs(:,1),ctrs(:,2),'ko',...
%     'MarkerSize',8,'LineWidth',2)
% guidata(hObject, handles); %updates the handles
% 
% 
% opts = statset('Display','final');
% [idx,ctrs] = kmeans(X67,2,...
%     'Distance','city',...
%     'Replicates',5,...
%     'Options',opts);
% axes(handles.axes14)
% scatter(X67(idx==1,1),X67(idx==1,2),'r');
% hold on
% scatter(X67(idx==2,1),X67(idx==2,2),'b');
% plot(ctrs(:,1),ctrs(:,2),'kx',...
%     'MarkerSize',8,'LineWidth',2)
% plot(ctrs(:,1),ctrs(:,2),'ko',...
%     'MarkerSize',8,'LineWidth',2)
% guidata(hObject, handles); %updates the handles
%**************************************************************************

        
%**************************************************************************         
% mean(F2_Normalization(dist56(2:end),'norm1'))
% mean(F2_Normalization(dist57(2:end),'norm1'))
% mean(F2_Normalization(dist67(2:end),'norm1'))
% 
% 'Pedram'
% 
% mean(dist56(2:end))
% mean(dist57(2:end))
% mean(dist67(2:end))
%**************************************************************************
% axes(handles.axes12)
% [z1, a1, b1, alpha1] = F0_fitellipse([pp(5,:);pp(6,:)], 'linear', 'constraint', 'trace');
% F0_plotellipse(z1, a1, b1, alpha1)
% 
% axes(handles.axes13)
% [z1, a1, b1, alpha1] = F0_fitellipse([pp(5,:);pp(7,:)],'linear', 'constraint', 'trace');
% F0_plotellipse(z1, a1, b1, alpha1)
% 
% 
% axes(handles.axes14)
% [z1, a1, b1, alpha1] = F0_fitellipse([pp(6,:);pp(7,:)], 'linear', 'constraint', 'trace');
% F0_plotellipse(z1, a1, b1, alpha1)

% [zl,rl]=F0_fitcircle([pp(6,:);pp(7,:)])
% t = linspace(0, 2*pi, 100);
% plot(zl(1),zl(2),'ro',zl(1) + rl * cos(t), zl(2)  + rl * sin(t), 'r')


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
handles.YY=YY;

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
%ylim([mean(HR1(1:tend))-50,mean(HR1(1:tend))+50]);
axis tight
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
%ylim([mean(BP1(1:tend))-50,mean(BP1(1:tend))+50]);
axis tight
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
        set(handles.popupmenu1,'String',{'AC02S','AA02S','AB02S','BH02S','CC02S','CG02S','DW02S','GB02S','GS02S','JF02S','JH02S','MSE02S','MSG02S','RM02S','RS02S','SB02S','SP02S','TH02S','FL02S','TC02S'});
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
cla(handles.axes11,'reset');
cla(handles.axes12,'reset');
cla(handles.axes13,'reset');
cla(handles.axes14,'reset');
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

axes(handles.axes11)
xlabel('Time [sec]','fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes12)
xlabel(param.title(5),'fontsize',10);
ylabel(param.title(6),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes13)
xlabel(param.title(5),'fontsize',10);
ylabel(param.title(7),'fontsize',10);
grid on
guidata(hObject, handles);

axes(handles.axes14)
xlabel(param.title(6),'fontsize',10);
ylabel(param.title(7),'fontsize',10);
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
        case 18
            file_name={'TH02S'};
        case 19
            file_name={'FL02S'};
        case 20
            file_name={'TC02S'};
            
    end
    rawdata=importdata(['../../../[Education] PhD/Dataset/ICORD/data/Rugby_games_HRV/dat02_cleaned/',cell2mat(file_name),'.dat']);
    x=rawdata.data;
    dl=length(x);
    if strcmp(file_name,'MSE02S')==1 || strcmp(file_name,'RM02S')==1 || strcmp(file_name,'CC02S')==1 || strcmp(file_name,'GS02S') || strcmp(file_name,'TC02S')==1
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
    plotATM1=['../../../[Education] PhD/Dataset/MIMIC/Data-MAT/',cell2mat(file_name),'nm.mat'];
    plotATM2=['../../../[Education] PhD/Dataset/MIMIC/Data-Raw/',cell2mat(file_name),'nm.info'];
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
CO(isnan(BP)) = 10;

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
if popup_sel_index12==2
    handles.data=[{HR},{BP},{Fs},{LoW},{dl},{CO}];
end
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
global param_mod;
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
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','c')
    xlabel(param_mod.title(1),'fontsize',10);
    ylabel('Frequency','fontsize',10);
    guidata(hObject, handles); %updates the handles
    grid on;
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes3,'reset');
    axes(handles.axes3)
    xlim([t1 t2]);
    hold on
    plot(t1:t1+size(pp,2)-1,pp(5,:),'c','LineWidth',2);
    grid on;
    box on;
    xlabel('Time [sec]','fontsize',10);
    ylabel(param_mod.title(1),'fontsize',10);
    guidata(hObject, handles); %updates the handles
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
global param;
global param_mod;
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
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','g')
    xlabel(param_mod.title(2),'fontsize',10);
    ylabel('Frequency','fontsize',10);
    guidata(hObject, handles); %updates the handles
    grid on;
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes5,'reset');
    axes(handles.axes5)
    xlim([t1 t2]);
    hold on
    plot(t1:t1+size(pp,2)-1,pp(6,:),'g','LineWidth',2);
    grid on;
    box on;
    xlabel('Time [sec]','fontsize',10);
    ylabel(param_mod.title(2),'fontsize',10);
    guidata(hObject, handles); %updates the handles
    
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
global param;
global param_mod;
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
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','b')
    xlabel(param_mod.title(3),'fontsize',10);
    ylabel('Frequency','fontsize',10);
    guidata(hObject, handles); %updates the handles
    grid on;
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes6,'reset');
    axes(handles.axes6)
    xlim([t1 t2]);
    hold on
    plot(t1:t1+size(pp,2)-1,pp(7,:),'b','LineWidth',2);
    grid on;
    box on;
    xlabel('Time [sec]','fontsize',10);
    ylabel(param_mod.title(3),'fontsize',10);
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


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global param;
LoW = param.LoW;
N_W = param.N_W;
YYtemp  = handles.YY;
YY=YYtemp';

YYmean= [mean(YYtemp,1);mean(YYtemp,1);mean(YYtemp,1)]'*1/3*ones(3,N_W);
YYmax = [max(YYtemp,[],1);max(YYtemp,[],1);max(YYtemp,[],1)]'*1/3*ones(3,N_W);
YYmin = [min(YYtemp,[],1);min(YYtemp,[],1);min(YYtemp,[],1)]'*1/3*ones(3,N_W);
YYnorm= (YY-YYmin)./(YYmax-YYmin);
% Hint: get(hObject,'Value') returns toggle state of checkbox8
if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    cla(handles.axes11,'reset');
    axes(handles.axes11)
    harea=area(YYnorm');
    set(harea(3),'FaceColor','c');
    grid on
    box on;
    guidata(hObject, handles); %updates the handles
    
else
    % Checkbox is not checked-take appropriate action
    cla(handles.axes11,'reset');
    axes(handles.axes11)
    harea=area(YYtemp);
    set(harea(3),'FaceColor','c');
    grid on;
    box on;
    guidata(hObject, handles); %updates the handles
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
