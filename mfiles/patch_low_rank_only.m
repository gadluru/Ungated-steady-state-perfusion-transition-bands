function update = patch_low_rank_only(Image,llr,para)
%--------------------------------------------------------------------------
%   update = patch_low_rank_only(Image,llr,para)
%--------------------------------------------------------------------------
%   computes the locally low rank update term for motion tracked patches
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,nset*nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames after cardiac phase fixing
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - llr: Motion tracked patches
%           - Npatch: number of motion tracked patches
%           - idx: linear indices of motion tracked patches [prod(patch_size),Ncycle,Npatch]
%           - mask_intensity: averaging mask of motion tracked patches [nx,ny,nset*nsl,nt]
%           - mask: location mask of motion tracked patches [nx,ny,nset*nsl,nt]
%       - para: reconstruction parameters [structure]
%           - Recon.tau: value for soft-thresholding
%--------------------------------------------------------------------------
%   Outputs:
%       - update: LLR update term [nx,ny,nt,nset*nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames after cardiac phase fixing
%           - nc: number of PCA coils
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%--------------------------------------------------------------------------

tau = para.Recon.tau;

Image = permute(Image,[1,2,4,3]); patch_all = Image(llr.idx);

patch_update = zeros(size(patch_all),'like',patch_all);
parfor i=1:llr.Npatch
    [U,S,V] = svd(patch_all(:,:,i),'econ');
    S = S - tau; S(S < 0) = 0;
    patch_update(:,:,i) = U*S*V';
end

update = accumarray(vec(llr.idx),vec(patch_update),[numel(Image),1]);
update = permute(reshape(update,size(Image)).*llr.mask_intensity,[1,2,4,3]); 

