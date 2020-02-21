function [] = fatigueX(n,varinow,DeltaN,inputs)
%defines inputs for strength model, and runs bundleX
%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global Xmax DX
global Lin Xavg CoV m Xin
global Lref Xref Lout
global Tsl T freebounds
global Df Vf
global k 
global N
global sigmath

%%%Input variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fatigue Cycles Defenition
%Variable jump 0cycle 

% DeltaN=zeros(100,1);
% DeltaN(1:100)=1;
% DeltaN(10001:15000)=500;
% DeltaN(2001:3000)=100;
% DeltaN(3001:3500)=1000; 
% DeltaN(3501:4000)=10000; 

%--------------------

N=length(DeltaN); % Number of iterations
R=inputs.R; % Fatigue Stress Ratio

%Numerical variables
Xmax=30000; DX=20;

%Input fibre strength distribution

Lin=inputs.Lin; Xavg=inputs.Xavg; CoV=inputs.CoV; 

%Interfacial shear strengthnn
Tsl=inputs.Tsl;

%Geometry and composite
Df=inputs.Df; Vf=inputs.Vf;
T=inputs.T; freebounds=inputs.freebounds;

Lref=inputs.Lref; 
% Output length
Lout=inputs.Lout; 

%Stress concentrations
k=varinow;

%%%Preliminary calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=inputs.m;   
% m=fzero(@(m) sqrt(gamma(1+2/m)/(gamma(1+1/m)^2)-1)-CoV, 1.2/CoV)
Xin=Xavg/(gamma(1+1/m));
Xref=Xin*((Lref/Lin)^(-1/m)); 

%%%Calculates strength distributions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bundleX(n,DeltaN,N,Df,R,inputs);
end