%% Fatigue code for unidirectional composites. 
% Developed by Marco Alves and Soraia Pimenta 2018
% Imperial College London,
% FiBreMoD Project 

%  More info can be found in the following publication:
% "A computationally-efficient micromechanical model for the fatigue life 
% of unidirectional composites under tension-tension loading"
% https://doi.org/10.1016/j.ijfatigue.2018.05.017


%%
%%clc
clear
clear all
close all
%% ============================= Inputs ===================================
% ======= fibre properties =======
inputs.E=240e3;%[MPa] 
inputs.Df=0.00685;%[mm] 
inputs.Vf=0.6; 
% fibre strength distribution
inputs.Lin=10;
inputs.Xavg=4293; 
inputs.m = 4.8;
inputs.CoV=0.25; 
% Reference length (unitary)
inputs.Lref=1;

%======= Specimen geometric properties =======
inputs.Lout=127; 
n=18; % bundle level ---> 2^i fibres

% ======= Matrix/Interface properties =======
% Interface shear strength
inputs.Tsl = 100; % MPa
%Toughness 
inputs.GIIc=2;%[kJ/m^2];
inputs.GIIthresh=0.0175;%[kJ/m^2];
% Paris law constants
inputs.C1=7.5e-5;
inputs.mp=2;

% ======= Shear lag boundary =======
%4-quadrangular, 6-hexagonal;
%1-interface, 3-matrix, 5-shortest.
inputs.T=41; inputs.freebounds=0;

% ========== Fatigue & numerical paramters ==========
inputs.R=0.1; % Fatigue Stress Ratio
nX = 1501; % number of stress increments - optimal value for computational efficiency

%Definition of the adaptive cycle step jump
DeltaCycles = [1 10 1000 ]; % cycle jump values
CycleIntervals = [100 3000 1000]; % number of increments for each cycle jump value

%%=========================== End of Inputs ===============================
for i = 1:length(DeltaCycles)
    Delta{i} = DeltaCycles(i)*ones(1,CycleIntervals(i));
end
% Pre-alocation
XavgoF = zeros(1,n+1,1);
leFF = zeros(nX,n+1,1);
lcontrolF = zeros(nX,n+1,1);
deltaAF = zeros(nX,n+1,1);

%%

for i = 1:length(DeltaCycles)
    DeltaN = Delta{i};
    fatigueX(n,2,DeltaN,inputs)
    global Fuo Xavgo X leF lcontrol deltaA DX
    if i==1
        XavgoF(:,:,1:size(Delta{i},2)) = Xavgo;
        leFF(:,:,1:size(Delta{i},2)) = leF;
        lcontrolF(:,:,1:size(Delta{i},2)) = lcontrol;
        deltaAF(:,:,1:size(Delta{i},2)) = deltaA;
        CyclesF(1) = DeltaN(1);
        for k =2:size(Delta{i},2)
            CyclesF(k) =CyclesF(k-1)+DeltaN(1); 
        end
    elseif i>1 && i<length(DeltaCycles)
        nextIndex = size(Delta{i-1},2) + size(Delta{i},2) - (size(Delta{i-1},2)*Delta{i-1}(1))/Delta{i}(1);
        aux = (size(Delta{i-1},2)*Delta{i-1}(1))/Delta{i}(1);   
        
        
        XavgoF(:,:,size(Delta{i-1},2)+1:nextIndex) = Xavgo(:,:,aux+1:size(Delta{i},2));
        leFF(:,:,size(Delta{i-1},2)+1:nextIndex) = leF(:,:,aux+1:size(Delta{i},2));
        lcontrolF(:,:,size(Delta{i-1},2)+1:nextIndex) = lcontrol(:,:,aux+1:size(Delta{i},2));
        deltaAF(:,:,size(Delta{i-1},2)+1:nextIndex) = deltaA(:,:,aux+1:size(Delta{i},2));
                CyclesF(size(Delta{i-1},2)+1) = CyclesF(size(Delta{i-1},2))+DeltaN(1);
        for k =size(Delta{i-1},2)+2:nextIndex
            CyclesF(k)=CyclesF(k-1)+DeltaN(1);
        end
    elseif i==length(DeltaCycles)
        [val, idx] = max(CyclesF);
        nextIndex = nextIndex+1;
      
        temp = 0;
        for j = 2:length(DeltaCycles)-1
            temp = temp + size(Delta{j},2) - (size(Delta{j-1},2)*Delta{j-1}(1))/Delta{j}(1);
        end
        aux1 = temp + size(Delta{1},2) +  size(Delta{end},2)- val/Delta{end}(1);
        aux2 = val/Delta{end}(1);
        XavgoF(:,:,nextIndex:aux1) = Xavgo(:,:,aux2+1:size(Delta{end},2));
        leFF(:,:,nextIndex:aux1) = leF(:,:,aux2+1:size(Delta{end},2));
        lcontrolF(:,:,nextIndex:aux1) = lcontrol(:,:,aux2+1:size(Delta{end},2));
        deltaAF(:,:,nextIndex:aux1) = deltaA(:,:,aux2+1:size(Delta{end},2));        
        CyclesF(nextIndex) = CyclesF(nextIndex-1)+DeltaN(1);
        for k =nextIndex+1:aux1
            CyclesF(k)=CyclesF(k-1)+DeltaN(1);
        end
    end
    
end
toc
%% =====================================================
%  ======================= Plots =======================
%  =====================================================
%% ======== Static Strength distribution - size effect
figure()
x0=20;
y0=15;
width=8;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
set(0,'DefaultTextInterpreter', 'latex')
Xm(:) = XavgoF(1,:,1);
plot([0:n],Xm)
xlabel('level -[$i$]');
ylab = ylabel('Strength (MPa)');
set(0,'DefaultTextInterpreter', 'latex')
set(gca,'TickLabelInterpreter', 'latex')

%% S-N Curves
%  Select Bundle level to be plotted:
bundle = 16;

% SN Curve
figure()
x0=20;
y0=15;
width=8;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
set(0,'DefaultTextInterpreter', 'latex')
strength(:) = XavgoF(1,bundle,:);
semilogx(CyclesF,strength/133000*0.6*100)
xlabel('Number of cycles');
ylab = ylabel('Stress (MPa)');
set(0,'DefaultTextInterpreter', 'latex')
set(gca,'TickLabelInterpreter', 'latex')

% Stochastic SN curve 
Xstatic=Xavgo(1,bundle,1)/DX;
maxStress = round(Xstatic);
minStress = 1;
staticexperimental=3270;
probabilityMatrix(:,:) = Fuo(minStress:maxStress,bundle,:)*100;
x = [1 CyclesF(end)];
y = [0 round(Xstatic*DX)]/1000;

figure()
x0=20;
y0=15;
width=9;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
set(0,'DefaultTextInterpreter', 'latex')
image(x,y,probabilityMatrix,'CDataMapping','scaled')
set(gca,'YDir','normal')
c=colorbar('TickLabelInterpreter','latex','FontSize',11);
ylabel(c,'Failure probability (\%)','interpreter','latex','FontSize',11)
colormap jet
set(0,'DefaultTextInterpreter', 'latex')
set(gca,'TickLabelInterpreter', 'latex')
xlabel ('Number of cycles','fontsize',11,'interpreter','latex')
ylabel ('Peak stress (GPa)','fontsize',11,'interpreter','latex')
set(gca,'yaxislocation');




