clc;
%0.1.1 Generalized method of computing transfer function of an equivalent
%0.1.2 transformer RLC ladder network.
%0.1.3 This is a driver code. For testing the effectiveness of the
%0.1.4 algorithm in the tester_code_underdampedInput_xxx.m file, first generate all
%0.1.5 processing data using this code (choose a particular value of the number 
%0.1.6 of coils in the transformer (n) then run the program for all 
%0.1.7 values of the location/node of partial discharge (i) from 2 to n 
%0.1.8 (both inclusive).
%1.0 Defining experimentally controlled variables
%1.0 Defining experimentally controlled variables
    %1.1.1 Number of RLC ladders (physically speaking, coils in the TF)
    %1.1.2 n should be between 3 and 7 (both included)
    n=4;
    %1.2.1 Partial discharge in the form of a current source Is
    %1.2.2 (to be defined later) Is=?
    %1.3.1 Location of node of partial discharge development (i):
    %1.3.2 Right now only ALONG the winding
    %1.2.3 i should be between 2 and n (both included)
    i=3;
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
%6.0 Plot obtained TF responses, etc.
%Define Partial Discharge input as a time varying current function, here,
%an underdamped sinusoidal signal
dt=1e-3;
t = 0:dt:1;
impulse= t==0;
A=5;f=10;%10kHz as time is in ms
x=A*exp(-1e0*abs(t)).*(1*sin(2*pi*f*t)+0*cos(2*pi*f*t));
%Save the Partial Discharge inputs for later use (verification)
filename=['x_time_domain_underdampedInput.mat'];
save(filename,'x');
%Show and save the time domain input partial discharge.
figure('units','normalized','outerposition',[0 0 1 1])
plot(t,x,'g+','LineWidth',2.5);
title(['Time Domain Graph for Input Underdamped Current (1 ms)']);
axis([0 1 -5 5]);
xlabel('Time');
ylabel('Amplitude');
ax = gca;
ax.FontSize = 13;
legend({'y = s'},'Location','northwest')
filename=['x_time_domain_underdampedInput'];
saveas(gcf,filename,'png')
pause(0.5);
%computing output values in frequency domain: 
%yL(s) = Line Short Circuit Current
%and yN(s) = Neutral Short Circuit Current
%Note that lsim function convolves the sysL/sysN (Transfer Function)
%with input (x(t)) in time domain, giving time domain short circuit
%Line/Neutral currents.
yL=(fft(lsim(sysL,x,t)));
yN=(fft(lsim(sysN,x,t)));
%The aim of this code was to compute (for given values of n, i and time domain
%partial discharge input function), the two output function yL(s) and yN(s)
%values and store them. Why? Since I am not actually experimenting on a real
%transformer, I can use these output values (as test inputs) to test the
%ability of the algorithm to detect the location of partial discharge for a
%given value of n and two sets of output values yL(t) and yLN(t) 
%(equivalently yL(s) and yN(s) in frequency domain).
len=length(yL);
q=-(len-1)/2:(len-1)/2;
filenameL=['yL_' num2str(n) '_' num2str(i) '_underdampedInput.mat'];
save(filenameL,'yL');
filenameN=['yN_' num2str(n) '_' num2str(i) '_underdampedInput.mat'];
save(filenameN,'yN');
%Computing Transfer Function values in frequency domain
xLimp=(fft(lsim(sysL,impulse,t)));
xNimp=(fft(lsim(sysN,impulse,t)));
%Applying Frequency Domain Relation (FDR) equation for both pairs of
%transfer functions and outputs. Therefore, theoritically speaking, both
%xL(s) and xN(s) are the fourier transform (frequency domain values) of the 
%same input partial discharge signal x(t).
xL=yL./xLimp;
xN=yN./xNimp;
%Obviously, Fourier Transform values cannot be plotted in 2D plane against
%frequency axis as the values are complex numbers. So we compute the
%magnitudes and plot them against frequency axis. 
%In other words, we plot the fourier spectrum of xL(s) and xN(s) 
fxL=sqrt(xL.*conj(xL));
fxN=sqrt(xN.*conj(xN));
%Show and save the predicted time domain frequency plots (magnitudes) of 
%partial discharge inputs.
plot(q,fxN,'--b','LineWidth',6)
hold on 
plot(q,fxL,':r','LineWidth',6)
hold off
title(['Fourier spectrum plots obained from transfer functions corresponding to n= ' num2str(n) ' and i= ' num2str(i) ' for pulse input']);
axis([0 2000 0 15]);
xlabel('Frequency');
ylabel('Amplitude');
ax = gca;
ax.FontSize = 13;
legend({'y = fxN','y = fxL'},'Location','northwest')
filename=['fxN_vs_fxL_' num2str(n) '_' num2str(i) '_underdampedInput'];
saveas(gcf,filename,'png')
%Evaluation of correctness of Partial Discharge predictions
correlation=corrcoef(xN,xL)
rmse_mtxN_mtxL=mean(sqrt((xN-xL).*conj(xN-xL)))
%***********************************************************************    