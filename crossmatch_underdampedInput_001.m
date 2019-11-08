clc
%0.1.1 Generalized method of computing transfer function of an equivalent
%0.1.2 transformer RLC ladder network.
%1.0 Defining experimentally controlled variables
    %1.1 Number of RLC ladders (physically speaking, coils in the TF)
    npotential=4;
    nactual=4;
    %1.2.1 Partial discharge in the form of a current source Is
    %1.2.2 (to be defined later) Is=?
    %1.3.1 Location of node of partial discharge development (i):
    %1.3.2 Right now only ALONG the winding
    ipotential=4;
    iactual=3;
%***********************************************************************    
%2.0 Defining ladder nw constants:
Rs=1.33;Cs=0.6;Cg=0.933;Ls=0.4310;
M=[0.2392,0.1435,0.0947,0.0496,zeros(1,npotential-5)];
%***********************************************************************    
%3.0 Constructing matrices useful for state space representation:
%3.1 Inductance (L) Matrix
L=zeros(npotential,npotential);
for r= 1:npotential
    for c= 1:npotential
        if r == c 
            L(r,c)=Ls;
        else
            L(r,c)=M(abs(r-c));
        end
    end
end
%3.2 Resistance (R) Matrix
R=Rs*eye(npotential);
%3.3.1 (S,SL1,SL2,SN1,SN2) Matrices
%3.3.2 L derived from when Live (node) shorted (to GND)
%3.3.3 N derived from when Neutral (node) shorted (to GND)
S=[-eye(npotential),zeros(npotential,1)]+[zeros(npotential,1),eye(npotential)];
SL1=S;SL1(:,1)=[];
SL2=zeros(npotential,1);SL2(ipotential-1,:)=1;
SN1=S;SN1(:,end)=[];
SN2=zeros(npotential,1);SN2(ipotential,:)=1;
%3.4 Capacitance (C) Matrix with spinoffs (CL1,CN1)
tmp=eye(npotential+1);tmp(1,1)=1/2;tmp(end,end)=1/2;
C=(Cg+2*Cs)*tmp+(-Cs)*[zeros(npotential,1),eye(npotential);zeros(1,npotential+1)]+(-Cs)*[zeros(1,npotential+1);eye(npotential),zeros(npotential,1)];
CL1=C;CL1(1,:)=[];CL1(:,1)=[];
CN1=C;CN1(end,:)=[];CN1(:,end)=[];
%***********************************************************************    
%4.0 Define (A),(B),(C),(D) matrices (State Space Representation):
    %But first, define (M),(-G),(T1):
    ML=[L,zeros(npotential,npotential);zeros(npotential,npotential),CL1];
    GL=-[-R,-SL1;SL1.',zeros(npotential,npotential)];
    T1L=[zeros(npotential,1);SL2];
    MN=[L,zeros(npotential,npotential);zeros(npotential,npotential),CN1];
    GN=-[-R,-SN1;SN1.',zeros(npotential,npotential)];
    T1N=[zeros(npotential,1);SN2];
AL=-ML\GL;
BL=ML\T1L;
AN=-MN\GN;
BN=MN\T1N;
tmpL=zeros(1,2*npotential);tmpL(:,1)=1;
tmpN=zeros(1,2*npotential);tmpN(:,npotential)=1;
CL=-Cs*AL(npotential+1,:)+tmpL;
DL=-Cs*BL(npotential+1,:);
CN=Cs*AN(2*npotential,:)+tmpN;
DN=Cs*BN(2*npotential,:);
%*********************************************************************** 
%5.0 Convert State Space Representation (SSR) into Transfer Function (TF)
[TFLb,TFLa]=ss2tf(AL,BL,CL,DL);
sysL=tf(TFLb,TFLa)
[TFNb,TFNa]=ss2tf(AN,BN,CN,DN);
sysN=tf(TFNb,TFNa)
%***********************************************************************
%6.0 Plot obtained TF responses, etc.
dt=1e-3;
t = 0:dt:1;
impulse= t==0;
testFilenameL=['yL_' num2str(nactual) '_' num2str(iactual) '_underdampedInput.mat'];
yLload=load(testFilenameL);
yL=[yLload.yL];
testFilenameN=['yN_' num2str(nactual) '_' num2str(iactual) '_underdampedInput.mat'];
yNload=load(testFilenameN);
yN=[yNload.yN];
xLimp=(fft(lsim(sysL,impulse,t)));
xNimp=(fft(lsim(sysN,impulse,t)));
len=length(xLimp);         %to take the frequency axis of the harmonics.
q=-(len-1)/2:(len-1)/2;  %divide the frequency compone
xL=yL./xLimp;
xN=yN./xNimp;
fxL=sqrt(xL.*conj(xL));
fxN=sqrt(xN.*conj(xN));
figure;
subplot(1,2,1);
plot(q,fxN,'-b')
hold on 
plot(q,fxL,':or')
hold off
title({
    ['Fourier Spectrum Plots predicted for'], 
    ['n = ' num2str(npotential) ' and i = ' num2str(ipotential)], 
    ['Actual values being n = ' num2str(nactual) ' and i = ' num2str(iactual)]
    });
axis([0 2000 0 15]);
xlabel('Frequency');
ylabel('Amplitude');
ax = gca;
ax.FontSize = 13;
legend({'y = xN','y = xL'},'Location','northwest')
filename=['xN_vs_xL_actual_' num2str(nactual) '_' num2str(iactual) '_potential_' num2str(npotential) '_' num2str(ipotential) '_underdampedInput'];
saveas(gcf,filename,'png') 
corrcoef(xN,xL)
rmse_xN_xL=mean(sqrt((xN-xL).*conj(xN-xL)))

txL=(ifft(xL));
mtxL=real(txL);
txN=(ifft(xN));
mtxN=real(txN);
subplot(1,2,2);
p1=plot(t,mtxN,'--b','LineWidth',6)
hold on 
p2=plot(t,mtxL,':r','LineWidth',6)
hold on
loadFilename=['x_time_domain_underdampedInput.mat'];
xload=load(loadFilename);
x=[xload.x];
p3=plot(t,x,'go')
p1.Color(4) = 1.00;
p2.Color(4) = 1.00;
p3.Color(4)= 0.05;
title({
    ['Partial Discharge Magnitudes in Time Domain predicted for' ] 
    ['n = ' num2str(npotential) ' and i = ' num2str(ipotential) ] 
    ['Actual values being n = ' num2str(nactual) ' and i = ' num2str(iactual)]
    });
axis([0 1 -10 10]);
xlabel('Time');
ylabel('Amplitude');
ax = gca;
ax.FontSize = 13;
legend({'y = mtxN','y = mtxL', 'y = x(t)'},'Location','northwest')
filename=['mtxN_vs_mtxL_actual_' num2str(nactual) '_' num2str(iactual) '_potential_' num2str(npotential) '_' num2str(ipotential) '_underdampedInput'];
saveas(gcf,filename,'png') 
corrcoef(mtxN,mtxL)
rmse_mtxN_mtxL=mean(sqrt((mtxN-mtxL).*conj(mtxN-mtxL)))