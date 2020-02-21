function [ leF, deltaA,lintact,gammath,lcontrol,Rglobal,iinit,idebond,iprop,lecrit,vCrackProp] = lefatigue( tau0, gamma0, Sinf, dSL, T, tm,  GIIth, GIIctm, G, gammaN, C1, m, lambda,  DeltaN, N, R, Lout,C,Ef,level,rf,A )
%Calculation of the effective recovery length increase due to the fatigue
%cycles for a specific value of remote stress

global DX
%initialize variables------------------------------------------------------
d=zeros(1,N+1);
gamma=zeros(1,N+1);
lcz=zeros(1,N);
deltaA=zeros(1,N);
erate=zeros(1,N);
dAdN=zeros(1,N);
dSL=2*Sinf;
iinit=0;
idebond=0;
iprop=0;
vCrackProp=0;

%Calculate initial parameters---------------------------------------
l=xDsigma(dSL,lambda,T,tau0); % first effective recovery length
le=zeros(1,N);
lcontrol=zeros(1,N);
lintact=zeros(1,N);
le(1)=l;
lintact(1)=le(1);
tauD(:,1)=ftau(l,tau0,lambda); % Shear stress in the critical point at N=1
damage(:,1)=1-tauD(:,1)/tau0;
gammath=gammathreshold( gammaN, tau0, GIIth, GIIctm); %Calculates gamma corresponding to GIIth
Rglobal=R;
lecrit=xtau(0,lambda,tau0);

%1st cycle iterations
gamma(1)=fgamma(le(1),gamma0, tau0, G, lambda ); % Gamma in the critical point


DeltaGII = ((Sinf)^2-(1-R^2))*T/Ef;
Gratio=DeltaGII/(GIIctm*tm);
vCrackProp=1*C1/C*(Gratio)^m; % crack increment - Paris law
%
Ninitiation = (lecrit-le(1))/vCrackProp;
DeltaNinitiation = ceil(Ninitiation/DeltaN(1))+1;


if gamma(1) < gammath % no damage can occur due to fatigue
    le(2:end)=le(1);
    lintact(2:end)=le(1);
else
    %Cycle Loop------------------------------------------------------------
    for i=2:N
        if le(i-1) < lecrit
            DeltaGII = (Sinf^2*(1-R^2))*T/Ef;
            Gratio=DeltaGII/(GIIctm*tm);
            vCrackProp=1*C1/C*(Gratio)^m; % crack increment - Paris law
            le(i)=le(i-1)+vCrackProp*DeltaN(i);
            lintact(i)=le(i);
        else
            le(i)=le(i-1)+vCrackProp*DeltaN(i);
            lintact(i:end)=le(i-1);
            deltaA(i)=le(i)-lecrit;
            for k = i+1:N
                le(k)=le(k-1)+vCrackProp*DeltaN(k);
                deltaA(k)=deltaA(k-1)+vCrackProp*DeltaN(k);
%                 if le(k)>Lout
%                     break
%                 end
            end
            break
        end
    end
end

%Correct the values of deltaA and lintact since the sum of both can never
%be bigger than Lout
boundary=(Lout/2)-(lecrit);
for i=1:N
    if deltaA(i) > boundary
        le(i:end)=Lout/2;
        lintact(i)=lintact(i)-(deltaA(i)-boundary);
        if lintact(i) < 0
            lintact(i:end)=0;
            deltaA(i:end)=Lout/2;
            for o=i:N
                idebond=idebond+DeltaN(o);
            end
            break
        end
    end
end

%Multiply by 2 because of the symmetry considered in the analysis
leF=2*le;
deltaA=2*deltaA;
lintact=2*lintact;

%Definition of the control length
lcontrol=2*leF;

for i=1:N
    if lcontrol(i)>Lout
        lcontrol(i)=Lout;
    end
end

end

