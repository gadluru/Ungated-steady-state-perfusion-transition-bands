function [Image,para] = STCR_conjugate_gradient_nSMS_ADMM(Data,para)
%--------------------------------------------------------------------------
%   [Image,para] = LLR_patch_tracking_4D_ADMM(Data,para)
%--------------------------------------------------------------------------
%   ADMM Step 1 Reconstruction for ungated, free-breathing, cardiac phase-resolved
%   myocardial perfusion MRI using Continuous Radial Interleaved
%   simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP). 
%
%   A standard cost function it solves is the spatially and temporally
%   constrained reconstruction (STCR):
%
%   || Am - d ||_2^2 + lambda_t || TV_t m ||_1 
%                    + lambda_t sum_i || P_i m ||_1
%                    + rho || Phat - m - Y ||_2^2
%
%   "A"         sampling matrix includes sensitivity maps, Fourier 
%               transform, and undersampling mask
%   "m"         image to be reconstructed
%   "d"         measured k-space data
%   ||.||_2^2   l2 norm
%   ||.||_1     l1 norm
%   "lambda_t"  temporal constraint weight
%   TV_t        temporal total variation (TV) operator (finite difference)
%               sqrt( abs(m_t+1 - m_t)^2 + epsilon )
%   "epsilon"   small term to aviod singularity
%   "rho"       l2 norm constraint weight
%   sum_i       sum across all patches
%   P_i         patch extraction operator
%   "Phat"      low-rank images
%   "Y"         ADMM Lagrangian multiplier
%--------------------------------------------------------------------------
%   Inputs:      
%       - Data: Reconstruction variables after cardiac phase fixing [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,nset,nsl]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt,1,1,nset,nsl]
%           - SMS: CAIPI phase modulation matrix [1,1,1,1,nsl,1,nsl]
%           - first_guess: prior STCR/LLR Step 3 Image (following first iteration of ADMM Step 1) [nx,ny,nt,nset*nsl]
%           - first_est: prior STCR [nx,ny,nt,nset*nsl]
%           - Y: ADMM Lagrangian multiplier (following first iteration of ADMM Step 1) [nx,ny,nt,nset*nsl]
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
%       - Image: ADMM Step 1 Reconstructed Perfusion Images [nx,ny,nt,nset*nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames after cardiac phase
%           fixing
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Reference:
%       [1] Whole-heart, ungated, free-breathing, cardiac-phase-resolved 
%           myocardial perfusion MRI by using Continuous Radial Interleaved
%           simultaneous Multi-slice acquisitions at sPoiled steady-state 
%           (CRIMP). MRM, in press.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

disp('Performing iterative CG reconstruction...');
disp('Showing progress...')

ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
weight_slTV    = para.Recon.weight_slTV;
weight_l2      = para.Recon.weight_l2;
epsilon        = para.Recon.epsilon;
para.Recon.step_size = para.step_size;

if isfield(Data,'first_est')
    new_img_x = single(Data.first_est);
else
    new_img_x = single(Data.first_guess);
end

for j=1:para.Recon.no_comp
    k = Data.kSpace(:,:,:,j,:,:,:);
    kSpace(:,j) = k(Data.mask);
end
Data.kSpace = kSpace; clear kSpace k

if ifGPU
    Data.kSpace   = gpuArray(Data.kSpace);
    new_img_x     = gpuArray(new_img_x);
    Data.sens_map = gpuArray(Data.sens_map);
    epsilon       = gpuArray(epsilon);
    if isfield(Data,'ramp_filter')
        Data.ramp_filter = gpuArray(Data.ramp_filter);
    end
    if isfield(Data,'sms_filter')
        Data.sms_filter = gpuArray(Data.sms_filter);
    end
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'SliceNorm',[],'l2Norm',[],'totalCost',[]);

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
temporal = @(im) compute_tTV_yt(im,weight_tTV,epsilon);
patch    = @(im) patch_ttv(im,Data.llr.MT,weight_tTV,epsilon);
bin      = @(im) compute_tTV_bins_no_lr(im,weight_tTV,para,epsilon);
spatial  = @(im) compute_sTV_yt(im,weight_sTV,epsilon);
slice    = @(im) compute_slTV_yt(im,weight_slTV,epsilon);

total_time = 0;
for iter_no = 1:para.Recon.noi

    if mod(iter_no,5) == 1
        t1 = tic;
    end

%% fidelity term/temporal/spatial TV
    [update_term,fidelity_norm] = fidelity(sr_forward(new_img_x,Data.size)); update_term = sr_backward(update_term,Data.size);
    
    update_term = update_term + temporal(new_img_x)*0.5;
    update_patch = patch(new_img_x)*0.5; update_bins = bin(new_img_x)*0.5;

    update_bins(permute(Data.llr.MT.mask,[1,2,4,3])) = update_patch(permute(Data.llr.MT.mask,[1,2,4,3]));
    update_term = update_term + update_bins; clear update_bins update_patch
    
    update_term = update_term + spatial(new_img_x);
    update_term = update_term + slice(new_img_x);
    
    if isfield(Data, 'Y')
        update_term = update_term + (Data.first_guess - new_img_x - Data.Y)*weight_l2;
    end

%% conjugate gradient
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+epsilon);
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%% line search
    if isfield(Data, 'Y')
        para.Cost = Cost_STCR_step_3(new_img_x, Data.first_guess - Data.Y, fidelity_norm, Data.llr.MT, para);
    else
        para.Cost = Cost_STCR_step_2(new_img_x, fidelity_norm, Data.llr.MT, para);
    end
    step_size = line_search_ADMM(new_img_x, update_term_old, Data, para);
    para.Recon.step_size(iter_no+1) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;

%% stoping creteria
    if para.Recon.break && iter_no > 1
        if step_size < 1e-4
            break
        end
    end
    
    if mod(iter_no,5) == 0
        fprintf(['Iteration = ' num2str(iter_no) '...']);
        total_time = total_time + toc(t1);
        toc(t1)
    end
end

Image = squeeze(gather(new_img_x));
figure, plotCost(para.Cost); drawnow
fprintf(['Iterative reconstruction running time is ' num2str(total_time) 's' '\n'])
end