function [Image,para] = LLR_patch_tracking_4D_ADMM(Data,para)
%--------------------------------------------------------------------------
%   [Image,para] = LLR_patch_tracking_4D_ADMM(Data,para)
%--------------------------------------------------------------------------
%   ADMM Reconstruction for ungated, free-breathing, cardiac phase-resolved
%   myocardial perfusion MRI using Continuous Radial Interleaved
%   simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP). 
%--------------------------------------------------------------------------
%   Inputs:      
%       - Data: Reconstruction variables after cardiac phase fixing [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,nset,nsl]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt,1,1,nset,nsl]
%           - SMS: CAIPI phase modulation matrix [1,1,1,1,nsl,1,nsl]
%           - first_guess: preliminary STCR [nx,ny,nt,nset*nsl]
%           - sens_map: coil sensitivity maps [nx,ny,1,nc,nsl,nset]
%           - ramp_filter: radial sampling filter [nx,ny,nt,1,1,nset]
%           - sms_filter: SMS acquisition filter [nx,ny,nt,1,1,nset]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames after cardiac phase fixing
%               - nc: number of PCA coils
%               - nset: number of slice groups
%               - nsl: number of sms slices (SMS, multiband=3)
%           - llr.MT: Motion tracked patches
%               - Npatch: number of motion tracked patches
%               - idx: linear indices of motion tracked patches [prod(patch_size),Ncycle,Npatch]
%               - mask_intensity: averaging mask of motion tracked patches [nx,ny,nset*nsl,nt]
%               - mask: location mask of motion tracked patches [nx,ny,nset*nsl,nt]
%           - llr.NM: Non-moving patches
%               - Npatch: number of non-moving patches [scalar]
%               - idx: linear indices of non-moving patches [prod(patch_size),Ncycle,Npatch]
%               - mask_intensity: averaging mask of non-moving patches [nx,ny,nset*nsl,nt]
%               - mask: location mask of non-moving patches [nx,ny,nset*nsl,nt]
%                   - nx: number of measurements in spatial x-dimension
%                   - ny: number of measurements in spatial y-dimension
%                   - nt: number of time frames after cardiac phase fixing
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%                   - Ncycle: number of cardiac cycles in the image series
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Image: Final ADMM Reconstructed Perfusion Images [nx,ny,nt,nset*nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames after cardiac phase
%           fixing
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------

para.Recon.weight_tTV = para.Recon.scale*para.weight_tTV;
para.Recon.weight_sTV = para.Recon.scale*para.weight_sTV;
para.Recon.weight_slTV = para.Recon.scale*para.weight_slTV;
para.Recon.noi = 20;

[Image,para] = STCR_conjugate_gradient_nSMS_ADMM(Data,para);

Data.Y = zeros(size(Data.first_guess),'single');
for ADMM_iter = 1:2
    %% ADMM update step 2
    % motion tracked patches
    Image_mc_lr = patch_low_rank_only(Image + Data.Y, Data.llr.MT, para);
    Image_mc_lr = permute(Image_mc_lr,[1,2,4,3]);
    
    % not moving patches
    Image_lr = LowRank_patch_yt(Image + Data.Y, Data.llr.NM, para);
    Image_lr = permute(Image_lr,[1,2,4,3]);

    % combine two LR images
    Image_lr(Data.llr.MT.mask) = Image_mc_lr(Data.llr.MT.mask);
    Image_lr = permute(Image_lr,[1,2,4,3]);
    
    %% ADMM update step 3
    Data.first_guess = Image_lr; Data.first_est = Image;
    Data.Y = Data.Y + Image - Image_lr;
    
    %% ADMM update step 1
    [Image,para] = STCR_conjugate_gradient_nSMS_ADMM(Data,para);
    
end

end

