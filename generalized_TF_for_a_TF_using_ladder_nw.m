%0.1.1 Generalized method of computing transfer function of an equivalent
%0.1.2 transformer RLC ladder network.
%1.0 Defining experimentally controlled variables
    %1.1 Number of RLC ladders (physically speaking, coils in the TF)
    n=4;
    %1.2.1 Partial discharge in the form of a current source Is
    %1.2.2 (to be defined later) Is=?
    %1.3.1 Location of node of partial discharge development (i):
    %1.3.2 Right now only ALONG the winding
    i=3;
%***********************************************************************    
%2.0 Defining ladder nw constants:
Rs=1.33;Cs=0.6;Cg=0.933;Ls=0.4310;
M=[0.2392,0.1435,0.0947,0.0612,0.0496,zeros(1,n-6)];
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
AL=-inv(ML)*GL;
BL=inv(ML)*T1L;
AN=-inv(MN)*GN;
BN=inv(MN)*T1N;
tmpL=zeros(1,2*n);tmpL(:,1)=1;
tmpN=zeros(1,2*n);tmpN(:,n)=1;
CL=[-Cs*AL(n+1,:)+tmpL];
DL=[-Cs*BL(n+1,:)];
CN=[Cs*AN(2*n,:)+tmpN];
DN=[Cs*BN(2*n,:)];
%*********************************************************************** 
%5.0 Convert State Space Representation (SSR) into Transfer Function (TF)
[TFLb,TFLa]=ss2tf(AL,BL,CL,DL);
sysL=tf(TFLb,TFLa)
[TFNb,TFNa]=ss2tf(AN,BN,CN,DN);
sysN=tf(TFNb,TFNa)
%***********************************************************************
%6.0 Plot obtained TF responses, etc.
subplot(4,2,1);
step(sysL);
subplot(4,2,2);
impulse(sysL);
subplot(4,2,3);
bode(sysL);
subplot(4,2,4);
nyquist(sysL);
grid;

subplot(4,2,5);
step(sysN);
subplot(4,2,6);
impulse(sysN);
subplot(4,2,7);
bode(sysN);
subplot(4,2,8);
nyquist(sysN);
grid;
%***********************************************************************    