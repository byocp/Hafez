%% Name: PEDRAM ATAEE             -            UBC Student Number: 32120073
%**************************************************************************
%                Fsolve Solution for Fowler's Model
%**************************************************************************
function out= F1_Fsolve_Sol(input)
global param_myfun;
param_myfun=input;
BP_init   = input(8);
HR_init   = input(9);
x0 = [HR_init; BP_init];              % Make a starting guess at the solution
options =optimset('Display','iter');  % Option to display output
options = optimset(options,'Display' ,'off');
options = optimset(options,'OutputFcn' ,{@outfun});
[x,fval,exitflag,output] = fsolve(@myfun,x0,options); % Call optimizer
global outfsolve;
Hf=x(1);
Pf=x(2);
out=[{Pf},{Hf},{outfsolve}];
end

function F = myfun(x)
global param_myfun;
%C_a      = param_myfun(1);
R0_a     = param_myfun(2);
DeltaV   = param_myfun(3);
IHR      = param_myfun(4);
%tau      = param_myfun(5);
V_H      = param_myfun(6);
Betta_H  = param_myfun(7);
%P_init   = param_myfun(8);
%HR_init  = param_myfun(9);
Alpha    = param_myfun(10);
gamma    = param_myfun(11);
Delta_h  = param_myfun(12);
sig_sp   = param_myfun(13);
sig_Alpha= param_myfun(14);

F=[x(1) - ((-1/Delta_h)*( V_H  * 1./(1+exp(-sig_Alpha.*(x(2)-sig_sp))) ...
    - Betta_H * (1-1./(1+exp(-sig_Alpha.*(x(2)-sig_sp)))) -Delta_h * IHR));
    x(2) - ((1+Alpha*(1-1./(1+exp(-sig_Alpha.*(x(2)-sig_sp))))) * R0_a * DeltaV * x(1))];
end

function stop=outfun(x,optimValues,state)
global outfsolve;
stop = false;
switch state
    case 'init'
        hold on
    case 'iter'
        outfsolve(optimValues.iteration+1,:)=x;
        %            history.fval = [history.fval; optimValues.fval];
        %            history.x = [history.x; x];
        %            searchdir = [searchdir;optimValues.searchdirection'];
        %            plot(x(1),x(2),'o');
        %            text(x(1)+.15,x(2),num2str(optimValues.iteration));
    case 'done'
        hold off
    otherwise
end
end