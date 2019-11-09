clc
%0.1.1 Generalized method of computing transfer function of an equivalent
%0.1.2 transformer RLC ladder network.
%0.1.3 I am assuming that you have already run the driver program
%(generator_code_underdampedInput_xxx.m file) for desired value(s) of n.
%1.0 Defining experimentally controlled variables
    %1.1 Number of RLC ladders (physically speaking, coils in the TF)
    %We aim to test whether our 'hypothesized' or 'potential' location/node of partial
    %discharge (ipotential) is in fact, the true/actual location/node of
    %partial discharge (iactual)
    nactual=4;
    %I'm assuming that we have prior information of the number of coils
    %in the transformer (i.e. we need not guess the value of n)
    npotential=nactual;
    %1.2.1 Partial discharge in the form of a current source Is
    %1.2.2 (to be defined later) Is=?
    %1.3.1 Location of node of partial discharge development (i):
    %1.3.2 Right now only ALONG the winding
    ipotential=3;
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
%6.0 Plot obtained TF responses, save output values, save plots etc.
%6.1.1 Define Partial Discharge input as a time varying current function, here,
%6.1.2 a pulse signal
dt=1e-3;
t = 0:dt:1;
impulse= t==0;
%6.2.1 Loading stored output values in frequency domain: 
%6.2.2 yL(s) = Line Short Circuit Current
%6.2.3 and yN(s) = Neutral Short Circuit Current
testFilenameL=['yL_' num2str(nactual) '_' num2str(iactual) '_underdampedInput.mat'];
yLload=load(testFilenameL);
yL=[yLload.yL];
testFilenameN=['yN_' num2str(nactual) '_' num2str(iactual) '_underdampedInput.mat'];
yNload=load(testFilenameN);
yN=[yNload.yN];
%6.3.1 Computing Transfer Function values in frequency domain
%6.2.4 Note that lsim function convolves the sysL/sysN (Transfer Function)
%6.2.5 with input (x(t)) in time domain, giving time domain short circuit
%6.2.6 Line/Neutral currents.
xLimp=(fft(lsim(sysL,impulse,t)));
xNimp=(fft(lsim(sysN,impulse,t)));
%6.4.1 Applying Frequency Domain Relation (FDR) equation for both pairs of
%6.4.2 transfer functions and outputs. Therefore, theoritically speaking, both
%6.4.3 xL(s) and xN(s) are the fourier transform (frequency domain values) of the 
%6.4.4 same input partial discharge signal x(t).
xL=yL./xLimp;
xN=yN./xNimp;
fxL=sqrt(xL.*conj(xL));
fxN=sqrt(xN.*conj(xN));
%6.4.9 Show the predicted time domain frequency plots (magnitudes) of 
%6.4.10 partial discharge inputs.
len=length(xLimp); 
q=-(len-1)/2:(len-1)/2;
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1);
plot(q,fxN,'--b','LineWidth',6)
hold on 
plot(q,fxL,':r','LineWidth',6)
hold off
title({
    ['Fourier spectrum plots predicted using transfer functions corresponding to n = ' num2str(npotential) ' and i = ' num2str(ipotential)],
    ['Actual number of coils being n = ' num2str(nactual) ' and actual location of partial discharge being i = ' num2str(iactual)]
    });
axis([0 2000 0 15]);
xlabel('Frequency');
ylabel('Amplitude');
ax = gca;
ax.FontSize = 13;
legend({'fxN(s) = Partial discharge signal predicted from short circuit neutral current','fxL(s) = Partial discharge signal predicted from short circuit line current'},'Location','northwest')
filename=['xN_vs_xL_actual_' num2str(nactual) '_' num2str(iactual) '_potential_' num2str(npotential) '_' num2str(ipotential) '_underdampedInput'];
saveas(gcf,filename,'png') 
%6.5.1 Applying Time Domain Relation (TDR) equation for both pairs of
%6.5.2 transfer functions and outputs. Therefore, theoritically speaking, both
%6.5.3 txL(t) and txN(t) are equal to the
%6.5.4 same input partial discharge signal x(t)
%6.5.5 There is no guarantee that a real signal subjected to FFT and then
%6.5.6 to IFFT will be come out the same. It could be that some portion of the 
%6.5.7 resultant signal may have imaginary parts. So we'll ignore these imaginary
%6.5.8 parts as errors
txL=(ifft(xL))
mtxL=real(txL);
txN=(ifft(xN));
mtxN=real(txN);
subplot(2,1,2);
p1=plot(t,mtxN,'--b','LineWidth',9)
hold on 
p2=plot(t,mtxL,':r','LineWidth',9)
hold on
loadFilename=['x_time_domain_underdampedInput.mat'];
xload=load(loadFilename);
x=[xload.x];
p3=plot(t,x,'g+');
p1.Color(4) = 1.00;
p2.Color(4) = 1.00;
p3.Color(4)= 0.05;
title({
    ['Partial discharge magnitudes in time domain predicted using transfer functions corresponding to n = ' num2str(npotential) ' and i = ' num2str(ipotential) ] 
    ['Actual number of coils being n = ' num2str(nactual) ' and actual location of partial discharge being i = ' num2str(iactual)]
    });
axis([0 1 -10 15]);
xlabel('Time');
ylabel('Amplitude');
ax = gca;
ax.FontSize = 13;
legend({'mtxN(t) = Partial discharge signal predicted from short circuit neutral current','mtxL(t) = Partial discharge signal predicted from short circuit line current','x(t) = Actual input partial discharge signal'},'Location','north')
filename=['mtxN_vs_mtxL_actual_' num2str(nactual) '_' num2str(iactual) '_potential_' num2str(npotential) '_' num2str(ipotential) '_underdampedInput'];
saveas(gcf,filename,'png') 
%7.0 Evaluation of correctness of Partial Discharge predictions
correlation_xN_xL=corrcoef(xN,xL)
rmse_xN_xL=mean(sqrt((xN-xL).*conj(xN-xL)))
correlation_mtxN_mtxL=corrcoef(mtxN,mtxL)
rmse_mtxN_mtxL=mean(sqrt((mtxN-mtxL).*conj(mtxN-mtxL)))
rmse_x_mtxN=mean(sqrt((mtxN-x').*conj(mtxN-x')))
rmse_x_mtxL=mean(sqrt((x'-mtxL).*conj(x'-mtxL)))