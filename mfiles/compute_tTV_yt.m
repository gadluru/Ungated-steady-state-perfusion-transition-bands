function tTV_update = compute_tTV_yt(image,weight,epsilon)
%--------------------------------------------------------------------------
%   tTV_update = compute_tTV_yt(image,weight,epsilon)
%--------------------------------------------------------------------------
%   computes the temporal total variation update term for time frames
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,nsl]
%       - weight: temporal total variation regularization weight [scalar]
%       - epsilon: small value to prevent singularity [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - update_ttv: temporal total vaiaration update term [nx,ny,nt,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------

if weight~=0
    temp_a = diff(image,1,3);
    temp_b = temp_a./(sqrt(epsilon+(abs(temp_a).^2)));
    temp_c = diff(temp_b,1,3);
    tTV_update = weight .* cat(3,temp_b(:,:,1,:,:,:),temp_c,-temp_b(:,:,end,:,:,:));
else
    tTV_update = 0;
end

end