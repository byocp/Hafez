%% Name: PEDRAM ATAEE             -            UBC Student Number: 32120073
%**************************************************************************
%  Lyapunov Analysis for Stability through Eigenvalues at Operating Point
%**************************************************************************
function out= F1_Lyap_EigVal(input)

C_a       = input(1);
R0_a      = input(2);
DeltaV    = input(3);
% IHR       = input(4);
tau       = input(5);
V_H       = input(6);
Betta_H   = input(7);
% P_init    = input(8);
% HR_init   = input(9);
Alpha     = input(10);
% gamma     = input(11);
Delta_h   = input(12);
sig_sp    = input(13);
sig_Alpha = input(14);
Pf        = input(15);

fcnsigd=@(x,y,z,coefd) (-1)*(-1)*x*coefd*exp(-1*x*(z-y))/(1+exp(-1*x*(z-y)))^2;
fcnsig=@(x,y,z)1/(1+exp(-1*x*(z-y)));

impmat22= ((-1)* R0_a * C_a * (1+Alpha*(1-fcnsig(sig_Alpha,sig_sp,Pf))) - ...
    (-1)*Pf * (-1) * R0_a * C_a * Alpha * fcnsigd(sig_Alpha,sig_sp,Pf,-1)) / ...
    (R0_a * C_a * (1+Alpha*(1-fcnsig(sig_Alpha,sig_sp,Pf))))^2;

impmat23= (-1)*(-1)*Pf * R0_a * C_a * Alpha * (-1)* fcnsigd(sig_Alpha,sig_sp,Pf,1) / ...
    (R0_a * C_a * (1+Alpha*(1-fcnsig(sig_Alpha,sig_sp,Pf))))^2;

impmat=[-1*Delta_h, fcnsigd(sig_Alpha,sig_sp,Pf,-1)* (-1)*(Betta_H)+fcnsigd(sig_Alpha,sig_sp,Pf,1)* (-1)*(V_H),  ...
    (-1)*Betta_H * fcnsigd(sig_Alpha,sig_sp,Pf,1);...
    DeltaV/C_a, impmat22, impmat23;...
    0,4/tau,-2/tau ];

eigvec=eig(impmat);
rleigvec=real(eigvec);
imeigvec=imag(eigvec);
out=[{rleigvec},{imeigvec}];

end

