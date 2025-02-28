function [Q,t] = Bullseye(Q)
%--------------------------------------------------------------------------
%   [Q,t] = Bullseye(Q)
%--------------------------------------------------------------------------
%   This function displays the blood flow quantification results using a
%   Bullseye plot.
%--------------------------------------------------------------------------
%   Inputs:      
%       - Q: Quantification parameters [struct]
%           -TCM
%               -para_low: lower bound of 2CM model parameters [1,4]
%               -para_up: upper bound of 2CM model parameters [1,4]
%               -para_init: initial estimate of 2CM model parameters [1,4]
%               -bld_flow_map: blood flow parameter map (identical to ktrans) [nx,ny,nset*nsl]
%               -kep_map: kep parameter map [nx,ny,nset*nsl]
%               -vb_map: vb parameter map [nx,ny,nset*nsl]
%               -dt_map: dt parameter map [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%           -mask
%               -Seg_mask: AHA regional binary mask for each slice [cell]
%               -Mask_all: myocardial binary mask for each slice [nx,ny,nset*nsl]
%               -Image_show: Display image for Bullseye [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%--------------------------------------------------------------------------
%   Outputs:
%       - Q: Quantification parameters [struct]
%           -TCM
%               -para_low: lower bound of 2CM model parameters [1,4]
%               -para_up: upper bound of 2CM model parameters [1,4]
%               -para_init: initial estimate of 2CM model parameters [1,4]
%               -bld_flow_map: blood flow parameter map (identical to ktrans) [nx,ny,nset*nsl]
%               -kep_map: kep parameter map [nx,ny,nset*nsl]
%               -vb_map: vb parameter map [nx,ny,nset*nsl]
%               -dt_map: dt parameter map [nx,ny,nset*nsl]
%               -bld_flow: Mean/SD AHA regional blood flow values [6,6,2]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%           -mask
%               -Seg_mask: AHA regional binary mask for each slice [cell]
%               -Mask_all: myocardial binary mask for each slice [nx,ny,nset*nsl]
%               -Image_show: Display image for Bullseye [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%       - t: Bullseye plot figure handle
%--------------------------------------------------------------------------
%   Reference:
%       [1] Adrian (2025). Bullseye Plot.zip (https://www.mathworks.com/
%           matlabcentral/fileexchange/47454-bullseye-plot-zip),
%           MATLAB Central File Exchange. Retrieved February 26, 2025. 
%--------------------------------------------------------------------------
Image = Q.mask.Image_show; Mask_all = Q.mask.Mask_all;

for i=1:size(Image,3)
    centroid(i,:) = regionprops(Mask_all(:,:,i)).Centroid;
end
x = round(centroid(:,2)); y = round(centroid(:,1));

for i=1:size(Image,3)
    temp(:,:,i) = Image(x(i)-80:x(i)+79,y(i)-80:y(i)+79,i);
end
Image = temp; clear temp

[nx,ny,nsl] = size(Image);

Image = reshape(Image,[nx,ny,nsl/3,3]);
Image = permute(Image,[1,3,2,4]);
Image = reshape(Image,[nx*(nsl/3),ny*3]);
Image = Image/max(Image(:));
Image = repmat(Image,[1,1,3]);

bld_flow_map = Q.TCM.bld_flow_map;

for i=1:size(bld_flow_map,3)
    temp(:,:,i) = bld_flow_map(x(i)-80:x(i)+79,y(i)-80:y(i)+79,i);
end
bld_flow_map = temp; clear temp

for i=1:size(Mask_all,3)
    temp(:,:,i) = Mask_all(x(i)-80:x(i)+79,y(i)-80:y(i)+79,i);
end
Mask_all = temp; clear temp

bld_map = bld_flow_map.*Mask_all;

bld_map = reshape(bld_map,[nx,ny,nsl/3,3]);
bld_map = permute(bld_map,[1,3,2,4]);
bld_map = reshape(bld_map,[nx*(nsl/3),ny*3]);

figure('Position',[0,0,1000,600]);

t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile(1);

imagesc(bld_map)
colormap hot
hold on
bg = imagesc(Image);

mask = logical(bld_map);
set(bg,'AlphaData',~mask);

colormap hot
axis image
axis off
caxis([0,2])

% Create the AHA 17-segment bullseye

for i=1:nsl
    Nseg = size(Q.mask.Seg_mask{i},3);
    bld_flow_map_temp = Q.TCM.bld_flow_map(:,:,i);
    for j=1:Nseg
        mask_temp = Q.mask.Seg_mask{i}(:,:,j);
        bld_flow_temp = bld_flow_map_temp(mask_temp);
        Q.TCM.bld_flow(i,j,1) = mean(bld_flow_temp);
        Q.TCM.bld_flow(i,j,2) = std(bld_flow_temp);
    end
end
bld_flow = Q.TCM.bld_flow;

bld_flow_mean = [];
for i=1:nsl
    bld_flow_slice = bld_flow_map(:,:,i);
    mask_slice = Mask_all(:,:,i);
    bld_flow_mean = [bld_flow_mean;vec(bld_flow_slice(mask_slice))];
end
bld_flow_std = std(nonzeros(bld_flow_mean));
bld_flow_mean = mean(nonzeros(bld_flow_mean));

number_of_regions = zeros(nsl,1);
for i=1:nsl
    number_of_regions(i) = length(nonzeros(bld_flow(i,:,1)));
end

fields = zeros(nsl+1,4); fields(1,:) = [0, 0.5, 1, 0];
for i=2:nsl+1
    fields(i,1:2) = fields(i-1,1:2) + 0.5;
    if number_of_regions(i-1) == 4
        fields(i,3:4) = [4,0];
    else
        fields(i,3:4) = [6,0];
    end
end

ax2 = nexttile(2);
c = createBullseye(fields); set(c,'Color','k','LineWidth',2), axis image, axis off

order = [2,1,6,5,4,3];
bld_flow = bld_flow(:,order,:);

% Example 1 of filling the bullseye, vector by vector
for i=1:nsl
    if number_of_regions(i) == 4
        fillBullseye(nonzeros(bld_flow(i,:,1)),fields(i+1,1),fields(i+1,2),0,360);
    else
        fillBullseye(nonzeros(bld_flow(i,:,1)),fields(i+1,1),fields(i+1,2),60,360+60);
    end
end

uistack(c,'top'),
colormap hot,

cbh2 = colorbar;
cbh2.TickLabels = num2cell(0:0.5:max([0,2]));
cbh2.Ticks = linspace(0,max([0,2]),length(cbh2.TickLabels));
caxis([0,2])

for i=1:nsl

    r = mean(fields(i+1,1:2));

    if number_of_regions(i) == 4
        theta = pi/4:pi/2:2*pi;

        x = cos(theta).*r; y = sin(theta).*r;

        rotation_all = [-45,45,-45,45]; idx = [1,2,5,6];
    else
        theta = (pi/3:pi/3:2*pi) + pi/6;

        x = cos(theta).*r; y = sin(theta).*r;

        rotation_all = [0,60,-60,0,60,-60]; idx = 1:6;
    end

    for j=1:length(idx)
        str = [num2str(nonzeros(bld_flow(i,idx(j),1)),'%0.2f'),'\pm',num2str(nonzeros(bld_flow(i,idx(j),2)),'%0.2f')];
        text(x(j),y(j),str,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',rotation_all(j),'FontSize',8);
    end
end
set(gca,'FontSize',12)

title(['Mean Flow = ',num2str(bld_flow_mean,'%0.2f'),'\pm',num2str(bld_flow_std,'%0.2f')],'FontSize',12)

end
