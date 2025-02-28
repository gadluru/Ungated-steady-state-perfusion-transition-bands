function [rv1,lv1]=auto_roi_mbr(cinemri)
%--------------------------------------------------------------------------
%   [rv1,lv1]=auto_roi_mbr(cinemri)
%--------------------------------------------------------------------------
%   This function returns a small region in the RV/LV blood pool to
%   generate an AIF for producing model-based images for registration
%--------------------------------------------------------------------------
%   Inputs:      
%       - cinemri: Unregistered Images [nx,ny,nt]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - rv1: RV blood pool location [nx,ny]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%       - lv1: LV blood pool location [nx,ny]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%--------------------------------------------------------------------------

[RV,LV] = Quant.FindLVRV(cinemri);

rv1=zeros(size(cinemri(:,:,15)));
lv1=zeros(size(cinemri(:,:,15)));

for i=-3:1:3
    for j=-3:1:3
        rv1(RV(1)+i+1,RV(2)+j+1)=1;
    end
end

for i=-3:1:3
    for j=-3:1:3
        lv1(LV(1)+i+1,LV(2)+j+1)=1;
    end
end


   
