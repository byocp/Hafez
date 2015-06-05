close all;
clear all;
load('results.mat')

figure(1)
%     subplot(3,2,1)
comp=[1 2 4];
incomp=[5 7 8];
LLa=[]; LLb=[];
dist56=[];
class_data_comp=[];
pedram56=[];
pedram56c=[];
for j=1:length(comp)
    i=comp(j);
    a=Imp.P56{1,i};
    pedram56=[pedram56;a];
    pedram56c=[pedram56c;ones(length(a),1)];
    hold on;
    grid on
    figure(j)
    scatter(a(:,1),a(:,2),'b^')
    class_data_comp=[class_data_comp;[a(:,1),a(:,2),ones(length(a),1)]]
%     dist56(i)=sqrt( sum(( a(:,1) - mean(a(:,1)) ).^2+( a(:,2) - mean(a(:,2)) ).^2) );
%     dist56(i)=sqrt(  sum( std(a(:,1)) + std(a(:,2))) .^2  );
    dist56(i,1)= std(a(:,1));
    dist56(i,2)= std(a(:,2));

end
class_data_incomp=[];
for j=1:length(incomp)
    i=incomp(j);
    a=Imp.P56{1,i};
    pedram56=[pedram56;a];
    pedram56c=[pedram56c;2*ones(length(a),1)];
    hold on;
        figure(j)
    scatter(a(:,1),a(:,2),'rv')
    class_data_incomp=[class_data_incomp; [a(:,1),a(:,2),2.*ones(length(a),1)]]
%     dist56(i)=sqrt( sum(( a(:,1) - mean(a(:,1)) ).^2+( a(:,2) - mean(a(:,2)) ).^2) );
%     dist56(i)=sqrt(  sum( std(a(:,1)) + std(a(:,2))) .^2  );
    dist56(i,1)= std(a(:,1));
    dist56(i,2)= std(a(:,2));
end

smplrnd1=randperm(length(class_data_comp));
smplrnd2=randperm(length(class_data_incomp));

class_data_train=[class_data_comp(smplrnd1( 1:floor(0.8*length(smplrnd1)) ),:); ...
    class_data_incomp(smplrnd2( 1:floor(0.8*length(smplrnd2)) ),:)];

class_data_test=[class_data_comp(smplrnd1( 1+floor(0.8*length(smplrnd1)):end ),:); ...
    class_data_incomp(smplrnd2( 1+floor(0.8*length(smplrnd2)):end ),:)];



X=class_data_test(:,1);
Y=class_data_test(:,2);
SL=class_data_train(:,1);
SW=class_data_train(:,2);
group=class_data_train(:,3);

[C,err,P,logp,coeff] = classify([X Y],[SL SW],...
    group,'quadratic');

figure(20)
% subplot(3,2,2)
grid on
% h1 = gscatter(SL,SW,group,'br','v^',[],'off');
% set(h1,'LineWidth',1)

hold on;
gscatter(X,Y,class_data_test(:,3),'br','*',3,'off');
legend('Complete','Incomplete','Location','NW')
K = coeff(1,2).const;
L = coeff(1,2).linear;
Q = coeff(1,2).quadratic;
f = sprintf('0 = %g+%g*x+%g*y+%g*x^2+%g*x.*y+%g*y.^2',...
    K,L,Q(1,1),Q(1,2)+Q(2,1),Q(2,2));
hold on
h2 = ezplot(f,[0 max(X) 0 max(Y)]);
set(h2,'Color','k','LineWidth',2)
title(['Classification with Test Data, Ac',num2str(err)])
xlabel('CVT');
ylabel('CST');
axis([0 max(X) 0 max(Y)])


figure(3)
gscatter(pedram56(:,1),pedram56(:,2),pedram56c,'br','^v');
hold on
legend('Complete','Incomplete')
% subplot(3,2,1)
h2 = ezplot(f,[0 max(X) 0 max(Y)]);
set(h2,'Color','k','LineWidth',2)
% title('Classification')
xlabel('Cardiac Parasympathetic Tone');
ylabel('Cardiac Sympathetic Tone');
title('');
axis([0 max(X) 0 max(Y)])
grid on



figure(1)
% subplot(3,2,3)
grid on
LLa=[]; LLb=[];
class_data_comp=[];
for j=1:length(comp)
    i=comp(j);
    a=Imp.P57{1,i};
    hold on;
    scatter(a(:,1),a(:,2),'b*')
    class_data_comp=[class_data_comp;[a(:,1),a(:,2),ones(length(a),1)]]
end
class_data_incomp=[];
for j=1:length(incomp)
    i=incomp(j);
    a=Imp.P57{1,i};
    hold on;
    scatter(a(:,1),a(:,2),'r*')
    class_data_incomp=[class_data_incomp; [a(:,1),a(:,2),2.*ones(length(a),1)]]
end

smplrnd1=randperm(length(class_data_comp));
smplrnd2=randperm(length(class_data_incomp));

class_data_train=[class_data_comp(smplrnd1( 1:floor(0.8*length(smplrnd1)) ),:); ...
    class_data_incomp(smplrnd2( 1:floor(0.8*length(smplrnd2)) ),:)];

class_data_test=[class_data_comp(smplrnd1( 1+floor(0.8*length(smplrnd1)):end ),:); ...
    class_data_incomp(smplrnd2( 1+floor(0.8*length(smplrnd2)):end ),:)];



X=class_data_test(:,1);
Y=class_data_test(:,2);
SL=class_data_train(:,1);
SW=class_data_train(:,2);
group=class_data_train(:,3);

[C,err,P,logp,coeff] = classify([X Y],[SL SW],...
    group,'quadratic');

figure(1)
subplot(3,2,4)
grid on
% h1 = gscatter(SL,SW,group,'br','v^',[],'off');
% set(h1,'LineWidth',1)
hold on;
gscatter(X,Y,class_data_test(:,3),'br','*',3,'off');
legend('Complete','Incomplete','Location','NW')
K = coeff(1,2).const;
L = coeff(1,2).linear;
Q = coeff(1,2).quadratic;
f = sprintf('0 = %g+%g*x+%g*y+%g*x^2+%g*x.*y+%g*y.^2',...
    K,L,Q(1,1),Q(1,2)+Q(2,1),Q(2,2));
hold on
h2 = ezplot(f,[0 max(X) 0 max(Y)]);
set(h2,'Color','k','LineWidth',2)

xlabel('')
ylabel('')
title('Classification with Test Data')
xlabel('CVT');
ylabel('AST');
axis([0 max(X) 0 max(Y)])

subplot(3,2,3)
h2 = ezplot(f,[0 max(X) 0 max(Y)]);
set(h2,'Color','k','LineWidth',2)
title('Classification')
xlabel('CVT');
ylabel('AST');
axis([0 max(X) 0 max(Y)])




figure(1)
subplot(3,2,5)
grid on
LLa=[]; LLb=[];
class_data_comp=[];
for j=1:length(comp)
    i=comp(j);
    a=Imp.P67{1,i};
    hold on;
    scatter(a(:,1),a(:,2),'b*')
    class_data_comp=[class_data_comp;[a(:,1),a(:,2),ones(length(a),1)]]
end
class_data_incomp=[];
for j=1:length(incomp)
    i=incomp(j);
    a=Imp.P67{1,i};
    hold on;
    scatter(a(:,1),a(:,2),'r*')
    class_data_incomp=[class_data_incomp; [a(:,1),a(:,2),2.*ones(length(a),1)]]
end

smplrnd1=randperm(length(class_data_comp));
smplrnd2=randperm(length(class_data_incomp));

class_data_train=[class_data_comp(smplrnd1( 1:floor(0.8*length(smplrnd1)) ),:); ...
    class_data_incomp(smplrnd2( 1:floor(0.8*length(smplrnd2)) ),:)];

class_data_test=[class_data_comp(smplrnd1( 1+floor(0.8*length(smplrnd1)):end ),:); ...
    class_data_incomp(smplrnd2( 1+floor(0.8*length(smplrnd2)):end ),:)];



X=class_data_test(:,1);
Y=class_data_test(:,2);
SL=class_data_train(:,1);
SW=class_data_train(:,2);
group=class_data_train(:,3);

[C,err,P,logp,coeff] = classify([X Y],[SL SW],...
    group,'quadratic');

figure(1)
subplot(3,2,6)
grid on
% h1 = gscatter(SL,SW,group,'br','v^',[],'off');
% set(h1,'LineWidth',1)
hold on;
gscatter(X,Y,class_data_test(:,3),'br','*',3,'off');
legend('Complete','Incomplete','Location','NW')
K = coeff(1,2).const;
L = coeff(1,2).linear;
Q = coeff(1,2).quadratic;
f = sprintf('0 = %g+%g*x+%g*y+%g*x^2+%g*x.*y+%g*y.^2',...
    K,L,Q(1,1),Q(1,2)+Q(2,1),Q(2,2));
hold on
h2 = ezplot(f,[0 max(X) 0 max(Y)]);
set(h2,'Color','k','LineWidth',2)

xlabel('')
ylabel('')
title('Classification with Test Data')
xlabel('CST');
ylabel('AST');
axis([0 max(X) 0 max(Y)])

subplot(3,2,5)
h2 = ezplot(f,[0 max(X) 0 max(Y)]);
set(h2,'Color','k','LineWidth',2)
title('Classification')
xlabel('CST');
ylabel('AST');
axis([0 max(X) 0 max(Y)])





figure(6)
h1 = gscatter(SL,SW,group,'br','v^',[],'off');
set(h1,'LineWidth',2)
legend('Complete','Incomplete','Location','NW')

hold on;
gscatter(X,Y,C,'br','*',3,'off');
K = coeff(1,2).const;
L = coeff(1,2).linear;
Q = coeff(1,2).quadratic;
f = sprintf('0 = %g+%g*x+%g*y+%g*x^2+%g*x.*y+%g*y.^2',...
    K,L,Q(1,1),Q(1,2)+Q(2,1),Q(2,2));
hold on
h2 = ezplot(f,[0 3 0 2]);
set(h2,'Color','k','LineWidth',2)
axis square
xlabel('')
ylabel('')
title('{\bf Classification with Fisher Training Data}')
xlabel('CST');
ylabel('AST');


figure(4)
subplot(1,3,1)
a=[]; aa=[]; aaa=[]; b=[]; bb=[]; bbb=[]; La=[]; Lb=[];
for j=1:length(comp)
    i=comp(j);
    a=[a;Imp.CVT{1,i}];
    aa=[aa;mean(Imp.CVT{1,i})];
    aaa=[aaa;std(Imp.CVT{1,i})];
end
for j=1:length(incomp)
    i=incomp(j);
    b=[b;Imp.CVT{1,i}];
    bb=[bb;mean(Imp.CVT{1,i})];
    bbb=[bbb;std(Imp.CVT{1,i})];
end
for i=1:length(a)
    La(i,:)='C';
end
for i=1:length(b)
    Lb(i,:)='I';
end
Label=char([La;Lb]);
Data=[a;b];
boxplot(Data,Label)
title('Cardiac Parasympathetic Tone');
[h,p,ci] = ttest2(aa,bb,[],[],'unequal')
kstest2(aa,bb);
[p,h]=signtest(aaa-bbb)




subplot(1,3,2)
a=[]; aa=[]; aaa=[]; b=[]; bb=[]; bbb=[]; La=[]; Lb=[];
for j=1:length(comp)
    i=comp(j);
    a=[a;Imp.CST{1,i}];
    aa=[aa;mean(Imp.CST{1,i})];
    aaa=[aaa;std(Imp.CST{1,i})];
end
for j=1:length(incomp)
    i=incomp(j);
    b=[b;Imp.CST{1,i}];
    bb=[bb;mean(Imp.CST{1,i})];
    bbb=[bbb;std(Imp.CST{1,i})];
end
for i=1:length(a)
    La(i,:)='C';
end
for i=1:length(b)
    Lb(i,:)='I';
end
Label=char([La;Lb]);
Data=[a;b];
boxplot(Data,Label)
title('Cardiac Sympathetic Tone');
[h,p,ci] = ttest2(aa,bb,[],[],'unequal')
kstest2(aa,bb)

subplot(1,3,3)
a=[]; aa=[]; aaa=[]; b=[]; bb=[]; bbb=[]; La=[]; Lb=[];
for j=1:length(comp)
    i=comp(j);
    a=[a;Imp.AST{1,i}];
    aa=[aa;mean(Imp.AST{1,i})];
    aaa=[aaa;std(Imp.AST{1,i})];
end
for j=1:length(incomp)
    i=incomp(j);
    b=[b;Imp.AST{1,i}];
    bb=[bb;median(Imp.AST{1,i})];
    bbb=[bbb;std(Imp.AST{1,i})];
end
for i=1:length(a)
    La(i,:)='C';
end
for i=1:length(b)
    Lb(i,:)='I';
end
Label=char([La;Lb]);
Data=[a;b];
boxplot(Data,Label)
title('Arterial Sympathetic Tone');
[h,p,ci] = ttest2(aa,bb,[],[],'unequal')
kstest2(aa,bb)


