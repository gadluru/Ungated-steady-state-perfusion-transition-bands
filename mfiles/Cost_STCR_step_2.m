function [Cost,totalCost] = Cost_STCR_step_2(Image, fNorm, llr, para)
%--------------------------------------------------------------------------
%   [Cost,totalCost] = Cost_STCR_step_2(Image, fNorm, llr, para)
%--------------------------------------------------------------------------
%   computes the cost of the ADMM Step 2 reconstruction problem 
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,nset*nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames after cardiac phase fixing
%               - nset: number of slice groups
%               - nsl: number of sms slices (SMS, multiband=3)
%       - fNorm: fidelity norm used for conjugate gradient step sizes [nk, nc]
%               - nk: number of acquired k-space samples
%               - nc: number of PCA coils
%       - llr: Motion tracked patches
%               - Npatch: number of motion tracked patches
%               - idx: linear indices of motion tracked patches [prod(patch_size),Ncycle,Npatch]
%               - mask_intensity: averaging mask of motion tracked patches [nx,ny,nset*nsl,nt]
%               - mask: location mask of motion tracked patches [nx,ny,nset*nsl,nt]
%                   - nx: number of measurements in spatial x-dimension
%                   - ny: number of measurements in spatial y-dimension
%                   - nt: number of time frames after cardiac phase fixing
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%                   - Ncycle: number of cardiac cycles in the image series
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Cost: variable containing cost of each regularization term 
%               for each iteration [structure]
%       - totalCost: total cost of the current step of the ADMM Step 2
%               reconstruction
%--------------------------------------------------------------------------

Cost = para.Cost; N = numel(Image);

fNorm = sum((abs(vec(fNorm)).^2)./prod(para.Recon.kSpace_size)/8);

if para.Recon.weight_tTV ~= 0
    tNorm1 = 0.5 .* para.Recon.weight_tTV .* sum(abs(vec(diff(Image,1,3))));

    Image_llr = permute(Image,[1,2,4,3]); patch_all = Image_llr(llr.idx); clear Image_llr
    tNorm2 = 0.5 .* para.Recon.weight_tTV .* sum(abs(vec(diff(patch_all,1,2))));

    tNorm = tNorm1 + tNorm2;
else
    tNorm = 0;
end

if para.Recon.weight_sTV ~= 0
    sx_norm = diff(Image,1,2); sx_norm(:,end+1,:,:,:) = 0;
    sy_norm = diff(Image,1,1); sy_norm(end+1,:,:,:,:) = 0;
    
    sNorm = para.Recon.weight_sTV .* sum(sqrt(vec(abs(sx_norm)).^2 + vec(abs(sy_norm)).^2));
else
    sNorm = 0;
end

if para.Recon.weight_slTV ~= 0
    slNorm = para.Recon.weight_slTV .* sum(abs(vec(diff(Image,1,4))));
else
    slNorm = 0;
end

slNorm = slNorm ./ N;
sNorm = sNorm ./ N;
tNorm = tNorm ./ N;
fNorm = fNorm ./ N;

totalCost = slNorm + sNorm + tNorm + fNorm;

if isempty(Cost.fidelityNorm)
    Cost.fidelityNorm = gather(fNorm);
    Cost.temporalNorm = gather(tNorm);
    Cost.spatialNorm = gather(sNorm);
    Cost.sliceNorm = gather(slNorm);
    Cost.totalCost = gather(totalCost);
else    
    Cost.fidelityNorm(end+1) = gather(fNorm);
    Cost.temporalNorm(end+1) = gather(tNorm);
    Cost.spatialNorm(end+1) = gather(sNorm);
    Cost.sliceNorm(end+1) = gather(slNorm);
    Cost.totalCost(end+1) = gather(totalCost);
end

end