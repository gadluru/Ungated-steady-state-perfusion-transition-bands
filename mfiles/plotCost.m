function plotCost(Cost)
%--------------------------------------------------------------------------
%   plotCost(Cost)
%--------------------------------------------------------------------------
%   plots the cost of each regularization term used for reconstruction for
%   each iteration
%--------------------------------------------------------------------------
%   Inputs:
%       - Cost: variable containing cost of each regularization term 
%               for each iteration [structure]
%--------------------------------------------------------------------------

loglog(Cost.totalCost,'LineWidth',2);
hold on;
plot(Cost.fidelityNorm,'LineWidth',2,'Marker','.','MarkerSize',15)
legend('Total Cost','Fidelity Norm')

if isfield(Cost,'temporalNorm')
    if ~isempty(Cost.temporalNorm)
        if Cost.temporalNorm(end)
            plot(Cost.temporalNorm,'x-','LineWidth',2,'MarkerSize',5);
            lgd = get(gca,'Legend');
            lgd = lgd.String;
            legend([lgd(1:end-1),'Temporal Norm'])
        end
    end
end

if isfield(Cost,'spatialNorm')
    if ~isempty(Cost.spatialNorm)
        if Cost.spatialNorm(end)
            plot(Cost.spatialNorm,'.-','LineWidth',2,'MarkerSize',15);
            lgd = get(gca,'Legend');
            lgd = lgd.String;
            legend([lgd(1:end-1),'Spatial Norm'])
        end
    end
end

if isfield(Cost,'sliceNorm')
    if ~isempty(Cost.sliceNorm)
        if Cost.sliceNorm(end)
            plot(Cost.sliceNorm,'o-','LineWidth',2,'MarkerSize',5);
            lgd = get(gca,'Legend');
            lgd = lgd.String;
            legend([lgd(1:end-1),'Slice Norm'])
        end
    end
end

if isfield(Cost, 'l2Norm')
    if ~isempty(Cost.l2Norm)
        if Cost.l2Norm(end)
            plot(Cost.l2Norm,'.-','LineWidth',2,'MarkerSize',15);
            lgd = get(gca,'Legend');
            lgd = lgd.String;
            legend([lgd(1:end-1),'L2 Norm'])
        end
    end
end

xlabel 'Iteration Number'
ylabel 'Norm'
set(gca,'FontSize',16)
hold off
