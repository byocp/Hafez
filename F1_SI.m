function out=F1_SI(vargin)
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