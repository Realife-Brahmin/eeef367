clc;
%0.1.1 Generalized method of computing transfer function of an equivalent
%0.1.2 transformer RLC ladder network.
%0.1.3 This is a driver code. For testing the effectiveness of the
%0.1.4 algorithm in the tester_code_pulseInput_xxx.m file, first generate all
%0.1.5 processing data using this code (choose a particular value of the number 
%0.1.6 of coils in the transformer (n) then run the program for all 
%0.1.7 values of the location/node of partial discharge (i) from 2 to n 
%0.1.8 (both inclusive).
%1.0 Defining experimentally controlled variables
    %1.1.1 Number of RLC ladders (physically speaking, coils in the TF)
    %1.1.2 n should be between 3 and 10 (both included)
    n=4;
    %1.2.1 Partial discharge is in the form of a current source Is
    %1.2.2 Location of node of partial discharge development (i):
    %1.2.3 Right now only ALONG the winding
    %1.2.4 i should be between 2 and n (both included)
    i=4;
%***********************************************************************    
%2.0 Defining ladder nw constants:
Rs=1.33;Cs=0.6;Cg=0.933;Ls=0.4310;
M=[0.2392,0.1435,0.0947,0.0496,zeros(1,n-5)];
%***********************************************************************    
%3.0 Constructing matrices useful for state space representation:
%3.1 Inductance (L) Matrix
L=zeros(n,n);
for r= 1:n
    for c= 1:n
        if r == c 
            L(r,c)=Ls;
        else
            L(r,c)=M(abs(r-c));
        end
    end
end
%3.2 Resistance (R) Matrix
R=Rs*eye(n);
%3.3.1 (S,SL1,SL2,SN1,SN2) Matrices
%3.3.2 L derived from when Live (node) shorted (to GND)
%3.3.3 N derived from when Neutral (node) shorted (to GND)
S=[-eye(n),zeros(n,1)]+[zeros(n,1),eye(n)];
SL1=S;SL1(:,1)=[];
SL2=zeros(n,1);SL2(i-1,:)=1;
SN1=S;SN1(:,end)=[];
SN2=zeros(n,1);SN2(i,:)=1;
%3.4 Capacitance (C) Matrix with spinoffs (CL1,CN1)
tmp=eye(n+1);tmp(1,1)=1/2;tmp(end,end)=1/2;
C=(Cg+2*Cs)*tmp+(-Cs)*[zeros(n,1),eye(n);zeros(1,n+1)]+(-Cs)*[zeros(1,n+1);eye(n),zeros(n,1)];
CL1=C;CL1(1,:)=[];CL1(:,1)=[];
CN1=C;CN1(end,:)=[];CN1(:,end)=[];
%***********************************************************************    
%4.0 Define (A),(B),(C),(D) matrices (State Space Representation):
    %But first, define (M),(-G),(T1):
    ML=[L,zeros(n,n);zeros(n,n),CL1];
    GL=-[-R,-SL1;SL1.',zeros(n,n)];
    T1L=[zeros(n,1);SL2];
    MN=[L,zeros(n,n);zeros(n,n),CN1];
    GN=-[-R,-SN1;SN1.',zeros(n,n)];
    T1N=[zeros(n,1);SN2];
AL=-ML\GL;
BL=ML\T1L;
AN=-MN\GN;
BN=MN\T1N;
tmpL=zeros(1,2*n);tmpL(:,1)=1;
tmpN=zeros(1,2*n);tmpN(:,n)=1;
CL=-Cs*AL(n+1,:)+tmpL;
DL=-Cs*BL(n+1,:);
CN=Cs*AN(2*n,:)+tmpN;
DN=Cs*BN(2*n,:);
%*********************************************************************** 
%5.0 Convert State Space Representation (SSR) into Transfer Function (TF)
[TFLb,TFLa]=ss2tf(AL,BL,CL,DL);
sysL=tf(TFLb,TFLa)
filenameL=['sysL_' num2str(n) '_' num2str(i)];
save(filenameL,'sysL');
[TFNb,TFNa]=ss2tf(AN,BN,CN,DN);
filenameN=['sysN_' num2str(n) '_' num2str(i)];
sysN=tf(TFNb,TFNa)
save(filenameN,'sysN');
%***********************************************************************
%6.0 Plot obtained TF responses, save output values, save plots etc.
%6.1.1 Define Partial Discharge input as a time varying current function, here,
%6.1.2 a pulse signal
dt=1e-3;
t = -1:dt:1;
impulse= t==0;
x=0.5*(sign(t+5*dt)-sign(t-5*dt));
%6.1.3 Save the Partial Discharge inputs for later use (verification)
filename=['x_time_domain_pulseInput.mat'];
save(filename,'x');
%6.1.4 Show and save the time domain input partial discharge plot.
figure('units','normalized','outerposition',[0 0 1 1])
plot(t,x,'g+','LineWidth',2.5);
title(['Time domain graph for input pulse current of width 10 ms']);
axis([-1e-1 1e-1 0 2]);
xlabel('Time');
ylabel('Amplitude');
ax = gca;
ax.FontSize = 13;
legend({'x(t) = Partial discharge current signal'},'Location','northwest')
filename=['x_time_domain_pulseInput'];
saveas(gcf,filename,'png')
pause(0.5);
%6.2.1 Computing output values in frequency domain: 
%6.2.2 yL(s) = Line Short Circuit Current
%6.2.3 and yN(s) = Neutral Short Circuit Current
%6.2.4 Note that lsim function convolves the sysL/sysN (Transfer Function)
%6.2.5 with input (x(t)) in time domain, giving time domain short circuit
%6.2.6 Line/Neutral currents.
yL=(fft(lsim(sysL,x,t)));
yN=(fft(lsim(sysN,x,t)));
%6.2.7 The aim of this code was to compute (for given values of n, i and time domain
%6.2.8 partial discharge input function), the two output function yL(s) and yN(s)
%6.2.9 values and store them. Why? Since I am not actually experimenting on a real
%6.2.10 transformer, I can use these output values (as test inputs) to test the
%6.2.11 ability of the walgorithm to detect the location of partial discharge for a
%6.2.12 given value of n and two sets of output values yL(t) and yLN(t) 
%6.2.13 (equivalently yL(s) and yN(s) in frequency domain).
filenameL=['yL_' num2str(n) '_' num2str(i) '_pulseInput.mat'];
save(filenameL,'yL');
filenameN=['yN_' num2str(n) '_' num2str(i) '_pulseInput.mat'];
save(filenameN,'yN');
%6.3.1 Computing Transfer Function values in frequency domain
xLimp=(fft(lsim(sysL,impulse,t)));
xNimp=(fft(lsim(sysN,impulse,t)));
%6.4.1 Applying Frequency Domain Relation (FDR) equation for both pairs of
%6.4.2 transfer functions and outputs. Therefore, theoritically speaking, both
%6.4.3 xL(s) and xN(s) are the fourier transform (frequency domain values) of the 
%6.4.4 same input partial discharge signal x(t).
xL=yL./xLimp;
xN=yN./xNimp;
%6.4.5 Obviously, Fourier Transform values cannot be plotted in 2D plane against
%6.4.6 frequency axis as the values are complex numbers. So we compute the
%6.4.7 magnitudes and plot them against frequency axis. 
%6.4.8 In other words, we plot the fourier spectrum of xL(s) and xN(s) 
fxL=sqrt(xL.*conj(xL));
fxN=sqrt(xN.*conj(xN));
%6.4.9 Show and save the predicted time domain frequency plots (magnitudes) of 
%6.4.10 partial discharge inputs.
len=length(yL);
q=-(len-1)/2:(len-1)/2;
plot(q,fxN,'--b','LineWidth',6)
hold on 
plot(q,fxL,':r','LineWidth',6)
hold off
title(['Fourier spectrum plots obained from transfer functions corresponding to n= ' num2str(n) ' and i= ' num2str(i) ' for pulse input']);
axis([0 2000 0 15]);%actual plot extends on both sides like the fft (sinc function) of a rect function should.
xlabel('Frequency');
ylabel('Amplitude');
ax = gca;
ax.FontSize = 13;
legend({'fxN(s) = Partial discharge signal predicted from short circuit neutral current','fxL(s) = Partial discharge signal predicted from short circuit line current'},'Location','northwest')
filename=['fxN_vs_fxL_' num2str(n) '_' num2str(i) '_pulseInput'];
saveas(gcf,filename,'png')
%***********************************************************************
%7.0 Evaluation of correctness of Partial Discharge predictions
correlation_xN_xL=corrcoef(xN,xL)
rmse_xN_xL=mean(sqrt((xN-xL).*conj(xN-xL)))
%***********************************************************************    