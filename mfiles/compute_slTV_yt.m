function slTV_update = compute_slTV_yt(Image,weight,epsilon)
%--------------------------------------------------------------------------
%   slTV_update = compute_slTV_yt(Image,weight,epsilon)
%--------------------------------------------------------------------------
%   computes the slice total variation update term
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,nsl]
%       - weight: slice total variation regularization weight [scalar]
%       - epsilon: small value to prevent singularity [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - slTV_update: slice total vaiaration update term [nx,ny,nt,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------

if weight~=0
    temp_a = diff(Image,1,4);
    temp_b = temp_a./(sqrt(epsilon+(abs(temp_a).^2)));
    temp_c = diff(temp_b,1,4);
    slTV_update = weight .* cat(4,temp_b(:,:,:,1),temp_c,-temp_b(:,:,:,end));
else
    slTV_update = 0;
end

end