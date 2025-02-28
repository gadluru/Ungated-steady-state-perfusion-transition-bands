
% This function produces the motion simulation results demonstrated in
% Supporting Information Figure 1 from the MRM publication,
% Mendes JK, Le JV, Arai AE, et al. Magn Reson Med 2025.
%--------------------------------------------------------------------------
%   Reference:
%       [1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G.
%           Improving ungated steady-state cardiac perfusion using transition
%           bands. Magn Reson Med 2025.
%--------------------------------------------------------------------------
%   Author:
%       Johnathan Le
%       E-mail: le.johnv@outlook.com
%--------------------------------------------------------------------------

addpath(genpath('../../mfiles/'))

T1 = 1:3000; % ms
TR = 2.4; % ms
flip_angle = 12; % degree
heart_rate_all = 120; % bpm
max_disp_all = 1; % slice thickness
Nalpha = 3000; %divisible by 6

Signal_TB_motion = zeros(6,Nalpha/6,length(T1),length(heart_rate_all),length(max_disp_all));
Signal_noTB_motion = zeros(6,Nalpha/6,length(T1),length(heart_rate_all),length(max_disp_all));
for i = 1:length(heart_rate_all)

    heart_rate = heart_rate_all(i);

    for j=1:length(max_disp_all)

        max_disp = max_disp_all(j);

        Signal_TB_motion(:,:,:,i,j) = get_TB_dic_with_phase_interleaving_2_SMS_3_with_motion(T1,flip_angle,TR,Nalpha,max_disp,heart_rate);
        Signal_noTB_motion(:,:,:,i,j) = get_noTB_dic_with_phase_interleaving_2_SMS_3_with_motion(T1,flip_angle,TR,Nalpha,max_disp,heart_rate);
    end
end

Signal_TB = get_TB_dic_with_phase_interleaving_2_SMS_3(T1,flip_angle,TR,Nalpha);
Signal_noTB = get_noTB_dic_with_phase_interleaving_2_SMS_3(T1,flip_angle,TR,Nalpha);

%%

TB_dic = squeeze(Signal_TB(1:3,end,:));
noTB_dic = squeeze(Signal_noTB(1:3,end,:));

Signal_TB = Signal_TB(1:3,:,:);
Signal_TB_motion = Signal_TB_motion(1:3,:,:);

Signal_noTB = Signal_noTB(1:3,:,:);
Signal_noTB_motion = Signal_noTB_motion(1:3,:,:);

SI_TB_tail = Signal_TB_motion(:,end/2:end,:);
SI_noTB_tail = Signal_noTB_motion(:,end/2:end,:);

SI_TB_max = squeeze(max(SI_TB_tail,[],2));
SI_TB_min = squeeze(min(SI_TB_tail,[],2));
SI_TB_mean = squeeze(mean(SI_TB_tail,2));

SI_noTB_max = squeeze(max(SI_noTB_tail,[],2));
SI_noTB_min = squeeze(min(SI_noTB_tail,[],2));
SI_noTB_mean = squeeze(mean(SI_noTB_tail,2));

for i=1:3
    d = abs(SI_noTB_max(i,:) - noTB_dic(i,:)');
    [~,maxCT1(i,:)] = min(d,[],1);

    d = abs(SI_noTB_min(i,:) - noTB_dic(i,:)');
    [~,minCT1(i,:)] = min(d,[],1);

    d = abs(SI_noTB_mean(i,:) - noTB_dic(i,:)');
    [~,meanCT1(i,:)] = min(d,[],1);

    d = abs(SI_TB_max(i,:) - TB_dic(i,:)');
    [~,maxFT1(i,:)] = min(d,[],1);

    d = abs(SI_TB_min(i,:) - TB_dic(i,:)');
    [~,minFT1(i,:)] = min(d,[],1);

    d = abs(SI_TB_mean(i,:) - TB_dic(i,:)');
    [~,meanFT1(i,:)] = min(d,[],1);
end

idx = 1500;

TB_motion_error = (abs(meanFT1(:,idx) - idx) ./ ((meanFT1(:,idx) + idx)/2)) * 100;
noTB_motion_error = (abs(meanCT1(:,idx) - idx) ./ ((meanCT1(:,idx) + idx)/2)) * 100;

fprintf('Transition Band Percent Motion Error - Slice 1: %.2f Slice 2: %.2f Slice 3: %.2f \n',TB_motion_error(1),TB_motion_error(2),TB_motion_error(3))
fprintf('No Transition Band Percent Motion Error - Slice 1: %.2f Slice 2: %.2f Slice 3: %.2f \n',noTB_motion_error(1),noTB_motion_error(2),noTB_motion_error(3))

%%

Colors = [0,0.447,0.7410; 0.85,0.325,0.098; 0.929,0.694,0.125; 0.494,0.184,0.556; 0.466,0.674,0.188; 0.301,0.745,0.933;0.6350,0.0780,0.1840];


total_time = (1:3000)*TR/1000; total_time = total_time(1:2:end); total_time = total_time(1:3:end);

fig = figure('Position',[0,0,1200,800]);

t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

%%

nexttile(4),plot(total_time,Signal_noTB(:,:,234)','LineWidth',1)
xlabel({'Time (s)';'\bf(d)'})
ylabel('Absolute SI [a.u.]'),ylim([60,120])
title('No Transition Band No Motion')
legend("Slice 1","Slice 2","Slice 3","Location","northwest")

nexttile(5),plot(total_time,Signal_noTB_motion(:,:,234)','LineWidth',1.5)
xlabel({'Time (s)';'\bf(e)'})
ylabel('Absolute SI [a.u.]'),ylim([60,120])
title('No Transition Band With Motion')
legend("Slice 1","Slice 2","Slice 3","Location","northwest")

nexttile(6),hold on

Y = [minCT1(1,:)',maxCT1(1,:)' - minCT1(1,:)'];
h1 = area(Y,'LineStyle','none');
h1(1).FaceAlpha = 0;
h1(2).FaceColor = (1-Colors(1,:))*0.6+Colors(1,:);

Y = [minCT1(2,:)',maxCT1(2,:)' - minCT1(2,:)'];
h2 = area(Y,'LineStyle','none');
h2(1).FaceAlpha = 0;
h2(2).FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:))*0.5+((1-Colors(1,:))*0.6+Colors(1,:))*0.5;
hlegend = area(3000,1,'LineStyle','none');
hlegend.FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:));

Y = [minCT1(3,:)',maxCT1(3,:)' - minCT1(3,:)'];
h3 = area(Y,'LineStyle','none');
h3(1).FaceAlpha = 0;
h3(2).FaceColor = (1-Colors(3,:))*0.6+Colors(3,:);

l0 = plot(T1,T1,'LineWidth',2,'Color','k');
l1 = plot(T1,meanCT1(1,:),'Color',Colors(1,:),'LineWidth',2);
l2 = plot(T1,meanCT1(2,:),'Color',Colors(2,:),'LineWidth',2);
l3 = plot(T1,meanCT1(3,:),'Color',Colors(3,:),'LineWidth',2);

legend([l0,l1,h1(2),l2,hlegend,l3,h3(2)],'Ideal','Slice 1','Slice 1 error','Slice 2','Slice 2 error','Slice 3','Slice 3 error','Location','northwest')
ylabel 'Fitted T1 (ms)'
xlabel({'Simulated T1 (ms)';'\bf(f)'})
title('No Transition Band T1 Motion Error')

axis([1,1500,1,1500])
xticks([1,300,600,900,1200,1500])
yticks([1,300,600,900,1200,1500,1800])

%%

nexttile(1),plot(total_time,Signal_TB(:,:,234)','LineWidth',1.5)
xlabel({'Time (s)';'\bf(a)'})
ylabel('Absolute SI [a.u.]'),ylim([60,120])
title('With Transition Band No Motion')
legend("Slice 1","Slice 2","Slice 3","Location","northwest")

nexttile(2),plot(total_time,Signal_TB_motion(:,:,234)','LineWidth',1.5)
xlabel({'Time (s)';'\bf(b)'})
ylabel('Absolute SI [a.u.]'),ylim([60,120])
title('With Transition Band With Motion')
legend("Slice 1","Slice 2","Slice 3","Location","northwest")

nexttile(3),hold on

Y = [minFT1(1,:)',maxFT1(1,:)' - minFT1(1,:)'];
h1 = area(Y,'LineStyle','none');
h1(1).FaceAlpha = 0;
h1(2).FaceColor = (1-Colors(1,:))*0.6+Colors(1,:);

Y = [minFT1(2,:)',maxFT1(2,:)' - minFT1(2,:)'];
h2 = area(Y,'LineStyle','none');
h2(1).FaceAlpha = 0;
h2(2).FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:))*0.5+((1-Colors(1,:))*0.6+Colors(1,:))*0.5;
hlegend = area(3000,1,'LineStyle','none');
hlegend.FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:));

Y = [minFT1(3,:)',maxFT1(3,:)' - minFT1(3,:)'];
h3 = area(Y,'LineStyle','none');
h3(1).FaceAlpha = 0;
h3(2).FaceColor = (1-Colors(3,:))*0.6+Colors(3,:);

l0 = plot(T1,T1,'LineWidth',2,'Color','k');
l1 = plot(T1,meanFT1(1,:),'Color',Colors(1,:),'LineWidth',2);
l2 = plot(T1,meanFT1(2,:),'Color',Colors(2,:),'LineWidth',2);
l3 = plot(T1,meanFT1(3,:),'Color',Colors(3,:),'LineWidth',2);

legend([l0,l1,h1(2),l2,hlegend,l3,h3(2)],'Ideal','Slice 1','Slice 1 error','Slice 2','Slice 2 error','Slice 3','Slice 3 error','Location','northwest')
ylabel 'Fitted T1 (ms)'
xlabel({'Simulated T1 (ms)';'\bf(c)'})
title('With Transition Band T1 Motion Error')

axis([1,1500,1,1500])
xticks([1,300,600,900,1200,1500])
yticks([1,300,600,900,1200,1500,1800])

exportgraphics(t,'MotionSimulation.png','Resolution',100)
