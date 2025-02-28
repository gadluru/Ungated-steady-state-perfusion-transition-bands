function tTV_update = compute_tTV_bins_no_lr(Image,weight,para,epsilon)
%--------------------------------------------------------------------------
%   tTV_update = compute_tTV_bins_no_lr(Image,weight,para,epsilon)
%--------------------------------------------------------------------------
%   computes the temporal total variation update term across cardiac cycles
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames after cardiac phase fixing
%           - nsl: total number of slices
%       - weight: temporal total variation regularization weight [scalar]
%       - para: reconstruction parameters [structure]
%           -Recon.bins: indices of different phases in the cardiac cycle [Nphase,nt]
%               -Nphase: Number of cardiac phases
%               -nt: number of time frames after cardiac phase fixing
%       - epsilon: small value to prevent singularity [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - update_ttv: temporal total vaiaration update term [nx,ny,nt,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames after cardiac phase fixing
%           - nsl: total number of slices
%--------------------------------------------------------------------------

tTV_update = zeros(size(Image),'like',Image); bins = para.Recon.bins;
for i=1:size(bins,1)
    tTV_update(:,:,bins(i,:),:) = compute_tTV_yt(Image(:,:,bins(i,:),:),weight,epsilon);
end