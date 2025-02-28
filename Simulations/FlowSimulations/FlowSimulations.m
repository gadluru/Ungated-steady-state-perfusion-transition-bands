
% This function produces the flow simulation results demonstrated in
% Figure 4 from the MRM publication, Mendes JK, Le JV, Arai AE, et al. 
% Magn Reson Med 2025.
%--------------------------------------------------------------------------
%   Reference:
%       [1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G.
%           Improving ungated steady-state cardiac perfusion using transition
%           bands. Magn Reson Med 2025.
%--------------------------------------------------------------------------
%   Author:
%       Jason Mendes
%       E-mail: Jason.Mendes@hsc.utah.edu
%--------------------------------------------------------------------------

SliceThickness=0.7; % cm
TR=2.39; % ms and every other beat
FlipAngle=12; % deg
PeakC=4; % mM
ShortC=1; % mM
load('PhantomData/Results.mat','noTB','TB_1','TB_3')
load('PhantomData/Velocity.mat','Velocity')

%% Tap water values

T1=2716;
Mz=1;
E1=exp(-2*TR./T1); % Every other pulse excitation
Mxy=zeros(1,1000);
for PulseIndex=1:1000
    Mxy(PulseIndex)=sin(12/180*pi)*Mz;
    Mz=Mz*cos(12/180*pi);
    Mz=1+(Mz-1).*E1;
end

% No Transition Band Simulation
[Ns,Nv,~]=size(noTB);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5),Ns,1,1),[1 Nv 1]);
SpinTravelPerExcite=repmat(Velocity*2*TR/1000,[Ns 1 1]);
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
noTB_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5+1),Ns,1,1),[1 Nv 1]);
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
TB_1_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5+2.4859),Ns,1,1),[1 Nv 1]); % Used transition of 2.4859 the slice thickness
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
TB_3_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);

% Correction for T1 data
f=noTB(:,:,1);
g=noTB_Sim;
f=min(f/max(f(:))*max(g(:)),max(Mxy));
Index_noTB=zeros(size(f));
for Index=1:numel(f)
    Index_noTB(Index)=find(Mxy<f(Index),1,'first');
end
PixelDiff=Mxy(Index_noTB-1)-Mxy(Index_noTB);
Weight_noTB=(Mxy(Index_noTB-1)-f)./PixelDiff;
f=TB_1(:,:,1);
g=TB_1_Sim;
f=f/max(f(:))*max(g(:));
Index_TB_1=zeros(size(f));
for Index=1:numel(f)
    Index_TB_1(Index)=find(Mxy<=f(Index),1,'first');
end
PixelDiff=Mxy(Index_TB_1-1)-Mxy(Index_TB_1);
Weight_TB_1=(Mxy(Index_TB_1-1)-f)./PixelDiff;
f=TB_3(:,:,1);
g=TB_3_Sim;
f=f/max(f(:))*max(g(:));
Index_TB_3=zeros(size(f));
for Index=1:numel(f)
    Index_TB_3(Index)=find(Mxy<=f(Index),1,'first');
end
PixelDiff=Mxy(Index_TB_3-1)-Mxy(Index_TB_3);
Weight_TB_3=(Mxy(Index_TB_3-1)-f)./PixelDiff;
noTB(:,:,1)=reshape(Mxy(Index_noTB).*Weight_noTB+Mxy(Index_noTB-1).*(1-Weight_noTB),Ns,Nv);
TB_1(:,:,1)=reshape(Mxy(Index_TB_1).*Weight_TB_1+Mxy(Index_TB_1-1).*(1-Weight_TB_1),Ns,Nv);
TB_3(:,:,1)=reshape(Mxy(Index_TB_3).*Weight_TB_3+Mxy(Index_TB_3-1).*(1-Weight_TB_3),Ns,Nv);

% noTB T1 map
M0=repmat(mean(noTB(:,1,1),1).*(1-E1*cos(12/180*pi))./(1-E1),[Ns Nv]);
dM=noTB(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
noTB_T1=-2*TR./log(E);
noTB_Error_T1=100*abs(T1-noTB_T1)./T1;
M0=repmat(noTB_Sim(:,1,1).*(1-E1*cos(12/180*pi))./(1-E1),[1 Nv]);
dM=noTB_Sim(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
noTB_Sim_T1=-2*TR./real(log(E));
noTB_Error_Sim_T1=100*abs(T1-noTB_Sim_T1)./T1;

%% Peak contrast values (4mM) r=3.7 prohance
fig1 = figure('Position',[0,0,1250,800]);
T1=63;
Mz=1;
E1=exp(-2*TR./T1); % Every other pulse excitation
Mxy=zeros(1,1000);
for PulseIndex=1:1000
    Mxy(PulseIndex)=sin(12/180*pi)*Mz;
    Mz=Mz*cos(12/180*pi);
    Mz=1+(Mz-1).*E1;
end

% noTB Simulation
[Ns,Nv,~]=size(noTB);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5),Ns,1,1),[1 Nv 1]);
SpinTravelPerExcite=repmat(Velocity*2*TR/1000,[Ns 1 1]);
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
noTB_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5+1),Ns,1,1),[1 Nv 1]);
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
TB_1_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5+2.4859),Ns,1,1),[1 Nv 1]); % Used transition of 2.4859 the slice thickness
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
TB_3_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);

% Correction for T1 data
noTB(:,:,1)=reshape(Mxy(Index_noTB).*Weight_noTB+Mxy(Index_noTB-1).*(1-Weight_noTB),Ns,Nv);
TB_1(:,:,1)=reshape(Mxy(Index_TB_1).*Weight_TB_1+Mxy(Index_TB_1-1).*(1-Weight_TB_1),Ns,Nv);
TB_3(:,:,1)=reshape(Mxy(Index_TB_3).*Weight_TB_3+Mxy(Index_TB_3-1).*(1-Weight_TB_3),Ns,Nv);

% noTB T1 map
M0=repmat(mean(noTB(:,1,1),1).*(1-E1*cos(12/180*pi))./(1-E1),[Ns Nv]);
dM=noTB(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
noTB_T1=-2*TR./log(E);
noTB_Error_T1=100*abs(T1-noTB_T1)./T1;
M0=repmat(noTB_Sim(:,1,1).*(1-E1*cos(12/180*pi))./(1-E1),[1 Nv]);
dM=noTB_Sim(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
noTB_Sim_T1=-2*TR./real(log(E));
noTB_Error_Sim_T1=100*abs(T1-noTB_Sim_T1)./T1;
subplot(2,3,1);
imagesc(noTB_Error_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Phantom)');
xlabel('Velocity (cm/s)');
title('No Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});
subplot(2,3,4);
imagesc(noTB_Error_Sim_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Simulated)');
xlabel('Velocity (cm/s)');
title('No Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});

% TB_1 T1 map
M0=repmat(mean(TB_1(:,1,1),1).*(1-E1*cos(12/180*pi))./(1-E1),[Ns Nv]);
dM=TB_1(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
TB_1_T1=-2*TR./log(E);
TB_1_Error_T1=100*abs(T1-TB_1_T1)./T1;
M0=repmat(TB_1_Sim(:,1,1).*(1-E1*cos(12/180*pi))./(1-E1),[1 Nv]);
dM=TB_1_Sim(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
TB_1_Sim_T1=-2*TR./real(log(E));
TB_1_Error_Sim_T1=100*abs(T1-TB_1_Sim_T1)./T1;
subplot(2,3,2);
imagesc(TB_1_Error_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Phantom)');
xlabel('Velocity (cm/s)');
title('7mm Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});
text(0.02,1.15,'PeakContrast (T1 = 63 ms)','Units','normalized','FontSize',12,'FontWeight','bold')
subplot(2,3,5);
imagesc(TB_1_Error_Sim_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Simulated)');
xlabel('Velocity (cm/s)');
title('7mm Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});

% TB_3 T1 map
M0=repmat(mean(TB_3(:,1,1),1).*(1-E1*cos(12/180*pi))./(1-E1),[Ns Nv]);
dM=TB_3(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
TB_3_T1=-2*TR./log(E);
TB_3_Error_T1=100*abs(T1-TB_3_T1)./T1;
M0=repmat(TB_3_Sim(:,1,1).*(1-E1*cos(12/180*pi))./(1-E1),[1 Nv]);
dM=TB_3_Sim(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
TB_3_Sim_T1=-2*TR./real(log(E));
TB_3_Error_Sim_T1=100*abs(T1-TB_3_Sim_T1)./T1;
subplot(2,3,3);
imagesc(TB_3_Error_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Phantom)');
xlabel('Velocity (cm/s)');
title('17.4mm Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});
subplot(2,3,6);
imagesc(TB_3_Error_Sim_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Simulated)');
xlabel('Velocity (cm/s)');
title('17.4mm Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});

%% Post injection values (1mM) r=3.7 prohance
fig2 = figure('Position',[0,0,1250,800]);
T1=234;
Mz=1;
E1=exp(-2*TR./T1); % Every other pulse excitation
Mxy=zeros(1,1000);
for PulseIndex=1:1000
    Mxy(PulseIndex)=sin(12/180*pi)*Mz;
    Mz=Mz*cos(12/180*pi);
    Mz=1+(Mz-1).*E1;
end

% noTB Simulation
[Ns,Nv,~]=size(noTB);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5),Ns,1,1),[1 Nv 1]);
SpinTravelPerExcite=repmat(Velocity*2*TR/1000,[Ns 1 1]);
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
noTB_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5+1),Ns,1,1),[1 Nv 1]);
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
TB_1_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);
SpinTravelDistance=repmat(reshape(SliceThickness*((Ns:-1:1)-0.5+2.4859),Ns,1,1),[1 Nv 1]); % Used transition of 2.4859 the slice thickness
ExcitePulses=min(SpinTravelDistance./(SpinTravelPerExcite+(SpinTravelPerExcite==0)*0.0001),length(Mxy));
Index=max(ceil(ExcitePulses),2);
Weight=1-(Index-ExcitePulses);
TB_3_Sim=reshape(Mxy(Index(:)).*Weight(:)'+Mxy(Index(:)-1).*(1-Weight(:)'),Ns,Nv);

% Correction for T1 data
noTB(:,:,1)=reshape(Mxy(Index_noTB).*Weight_noTB+Mxy(Index_noTB-1).*(1-Weight_noTB),Ns,Nv);
TB_1(:,:,1)=reshape(Mxy(Index_TB_1).*Weight_TB_1+Mxy(Index_TB_1-1).*(1-Weight_TB_1),Ns,Nv);
TB_3(:,:,1)=reshape(Mxy(Index_TB_3).*Weight_TB_3+Mxy(Index_TB_3-1).*(1-Weight_TB_3),Ns,Nv);

% noTB T1 map
M0=repmat(mean(noTB(:,1,1),1).*(1-E1*cos(12/180*pi))./(1-E1),[Ns Nv]);
dM=noTB(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
noTB_T1=-2*TR./log(E);
noTB_Error_T1=100*abs(T1-noTB_T1)./T1;
M0=repmat(noTB_Sim(:,1,1).*(1-E1*cos(12/180*pi))./(1-E1),[1 Nv]);
dM=noTB_Sim(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
noTB_Sim_T1=-2*TR./real(log(E));
noTB_Error_Sim_T1=100*abs(T1-noTB_Sim_T1)./T1;
subplot(2,3,1);
imagesc(noTB_Error_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Phantom)');
xlabel('Velocity (cm/s)');
title('No Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});
subplot(2,3,4);
imagesc(noTB_Error_Sim_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Simulated)');
xlabel('Velocity (cm/s)');
title('No Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});

% TB_1 T1 map
M0=repmat(mean(TB_1(:,1,1),1).*(1-E1*cos(12/180*pi))./(1-E1),[Ns Nv]);
dM=TB_1(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
TB_1_T1=-2*TR./log(E);
TB_1_Error_T1=100*abs(T1-TB_1_T1)./T1;
M0=repmat(TB_1_Sim(:,1,1).*(1-E1*cos(12/180*pi))./(1-E1),[1 Nv]);
dM=TB_1_Sim(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
TB_1_Sim_T1=-2*TR./real(log(E));
TB_1_Error_Sim_T1=100*abs(T1-TB_1_Sim_T1)./T1;
subplot(2,3,2);
imagesc(TB_1_Error_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Phantom)');
xlabel('Velocity (cm/s)');
title('7mm Transition Band');
text(0.02,1.15,'PostContrast (T1 = 234 ms)','Units','normalized','FontSize',12,'FontWeight','bold')
xticks(1:5:16);
xticklabels({'0','10','20','30'});
subplot(2,3,5);
imagesc(TB_1_Error_Sim_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Simulated)');
xlabel('Velocity (cm/s)');
title('7mm Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});

% TB_3 T1 map
M0=repmat(mean(TB_3(:,1,1),1).*(1-E1*cos(12/180*pi))./(1-E1),[Ns Nv]);
dM=TB_3(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
TB_3_T1=-2*TR./log(E);
TB_3_Error_T1=100*abs(T1-TB_3_T1)./T1;
M0=repmat(TB_3_Sim(:,1,1).*(1-E1*cos(12/180*pi))./(1-E1),[1 Nv]);
dM=TB_3_Sim(:,:,1)./M0;
E=(1-dM)./(1-dM.*cos(12/180*pi));
TB_3_Sim_T1=-2*TR./real(log(E));
TB_3_Error_Sim_T1=100*abs(T1-TB_3_Sim_T1)./T1;
subplot(2,3,3);
imagesc(TB_3_Error_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Phantom)');
xlabel('Velocity (cm/s)');
title('17.4mm Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});
subplot(2,3,6);
imagesc(TB_3_Error_Sim_T1(:,1:16),[0 100]);
h = colorbar;
set(get(h,'title'),'String','% Error');
ylabel('Slice Number (Simulated)');
xlabel('Velocity (cm/s)');
title('17.4mm Transition Band');
xticks(1:5:16);
xticklabels({'0','10','20','30'});

%%

exportgraphics(fig1,'PeakContrast_FlowSimulation.png','Resolution',100)
exportgraphics(fig2,'PostContrast_FlowSimulation.png','Resolution',100)

