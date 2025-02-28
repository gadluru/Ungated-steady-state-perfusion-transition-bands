function [RV, LV] = FindLVRV(cinemri1)
%--------------------------------------------------------------------------
%   [RV, LV] = FindLVRV(cinemri1)
%--------------------------------------------------------------------------
%   This function returns the pixel location of a region of interest in the
%   RV/LV blood pool to generate an AIF for producing model-based images
%   for registration.
%--------------------------------------------------------------------------
%   Inputs:      
%       - cinemri1: Unregistered Images [nx,ny,nt]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - RV: RV blood pool pixel location [1,2]
%       - LV: LV blood pool pixel location [1,2]
%--------------------------------------------------------------------------

[nx ny nt] = size(cinemri1);

smoothedcinemri1 = imgaussfilt(cinemri1,1);

mymax = squeeze(max(smoothedcinemri1,[],3));

smoothed = imgaussfilt(mymax,4);

BW = imregionalmax(smoothed);

[x, y] = find(BW);

temporalSmoothing = -diff(fspecial('gaussian',[1 (nt+10+1)],3));
maxUpslope = zeros(length(x),3);

for i=1:length(x)
    mask_temp = false(nx,ny);
    mask_temp(x(i),y(i)) = true;
    mask_temp = bwdist(mask_temp) < 5;
    N_pixel_in_mask = sum(sum(mask_temp));
    values(i,:) = sum(sum(cinemri1 .* mask_temp))/N_pixel_in_mask;
    mean_first_5 = mean(values(i,1:5));
    mean_last_5 = mean(values(i,end-4:end));
    values(i,:) = values(i,:) - mean_first_5;
    signal = [repmat(mean_first_5,[1,5]),values(i,:),repmat(mean_last_5,[1,5])];
    smoothedCurve = filter2(temporalSmoothing,signal);
    [maxUpslope(i,1) maxUpslope(i,2)] = max(smoothedCurve);
    maxUpslope(i,3) = norm([(x(i) - nx/2) (y(i) - ny/2)]); 
end

score = exp(-(maxUpslope(:,3)).^2/(2*(nx/12)^2)).*(1-1./(1+(maxUpslope(:,1)/(max(maxUpslope(:,1)))/2).^2));
if(length(score) >2)
    IDX = score>(max(score)/2);
    ratio = 3;
    while(sum(IDX) < 2)
        IDX = score>(max(score)/ratio);
        ratio = ratio + 1;
    end
    [~, worst] = min(score);
    NotWanted = (IDX == IDX(worst));

    x(NotWanted) = [];
    y(NotWanted) = [];

    maxUpslope(NotWanted,:) = [];
    score(NotWanted) = [];
end
IDX = kmeans(maxUpslope(:,2),2,'EmptyAction','singleton');

RVspot = find(IDX == 1);
weights = score(RVspot);
weights = weights / sum(weights);
tempx = x(RVspot);
tempy = y(RVspot);

RV = round([sum(tempx.*weights) sum(tempy.*weights)]);
LVspot = find(IDX == 2);
weights = score(LVspot);
weights = weights / sum(weights);
tempx = x(LVspot);
tempy = y(LVspot);

LV = round([sum(tempx.*weights) sum(tempy.*weights)]);

if(LV(2) < RV(2))
    temp = RV;
    RV = LV;
    LV = temp;
end

end