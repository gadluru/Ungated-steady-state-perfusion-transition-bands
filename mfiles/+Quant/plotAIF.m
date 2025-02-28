function plotAIF(Gd_curve,TimeStamps)
clrs = [0,0.447,0.7410; 0.85,0.325,0.098; 0.929,0.694,0.125; 0.494,0.184,0.556; 0.466,0.674,0.188; 0.301,0.745,0.933;0.6350,0.0780,0.1840];

labels = [];
for i=1:size(Gd_curve,2)
    labels = [labels,strcat("Slice",num2str(i))];
end

figure,hold on
for i=1:size(Gd_curve,2)
    plot(TimeStamps,Gd_curve(:,i),'LineWidth',2,'Color',clrs(i,:))
end
title ('Arterial Input Function')
xlabel 'Time [second]'
ylabel 'Gd concentration [mmol/L]'
legend(labels)