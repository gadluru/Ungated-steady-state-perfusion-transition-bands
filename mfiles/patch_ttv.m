function update = patch_ttv(Image,llr,weight,epsilon)
%--------------------------------------------------------------------------
%   update = patch_ttv(Image,llr,weight,epsilon)
%--------------------------------------------------------------------------
%   computes the patch-based temporal total variation update term
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,nsl]
%       - llr: Motion tracked patches
%           - Npatch: number of motion tracked patches
%           - idx: linear indices of motion tracked patches [prod(patch_size),Ncycle,Npatch]
%           - mask_intensity: averaging mask of motion tracked patches [nx,ny,nset*nsl,nt]
%           - mask: location mask of motion tracked patches [nx,ny,nset*nsl,nt]
%       - weight: temporal total variation regularization weight [scalar]
%       - epsilon: small value to prevent singularity [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - update: patch-based temporal total vaiaration update term [nx,ny,nt,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames after cardiac phase fixing
%           - nsl: total number of slices
%--------------------------------------------------------------------------

Image = permute(Image,[1,2,4,3]);
patch_all = Image(llr.idx);

update_tv = compute_tTV_yt(permute(patch_all,[1,3,2]),weight,epsilon);
update_tv = permute(update_tv,[1,3,2]);

update_tv = accumarray(vec(llr.idx),vec(update_tv),[numel(Image),1]);
update_tv = reshape(update_tv,size(Image));

update = update_tv.*llr.mask_intensity;    
update = permute(update,[1,2,4,3]);

end