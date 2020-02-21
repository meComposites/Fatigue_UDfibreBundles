function [ tau0, gamma0, dSL, GIIth, GIIctm, G, gammaN, C1, mp,  lambda,E ] = fproperties( X, Tf,tm,Tsl,inputs )
%Marco Alves

%==============Fiber Young Modulus==============
E=inputs.E;%[MPa] % T700 CF
%===============================================



%==============B.C @ xi=0==============
dS0=0; %[Pa];
tau0=Tsl;%[MPa];
gamma0=0;
%============BC @ xi/le/2=============
dSL=2*X;%[MPa];
%======================================

%==========Matrix Properties===========
GIIc=inputs.GIIc;%[kJ/m^2];
GIIthresh=inputs.GIIthresh;
GIIctm=GIIc/tm;
gammaN=(GIIctm)*2/(tau0);
GIIth=GIIthresh/tm;
G=-(tau0/gammaN);%[MPa];
%======================================

%==========Fatigue Parameters==========
% Paris law constants
C1=inputs.C1;
mp=inputs.mp;

% Shear lag paramters
x0=0;
lambda=sqrt(2*abs(G)./(Tf*tm*E));
%======================================
end

