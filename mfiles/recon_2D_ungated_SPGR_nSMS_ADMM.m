function recon_2D_ungated_SPGR_nSMS_ADMM(kSpace_all,para)
%--------------------------------------------------------------------------
%   recon_2D_ungated_SPGR_mSMS_ADMM(kSpace_all,para)
%--------------------------------------------------------------------------
%   Reconstruction for ungated, free-breathing, cardiac phase-resolved
%   myocardial perfusion MRI using Continuous Radial Interleaved
%   simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP). 
%   This reconstruction code builds upon the CRIMP framework
%   to perform reconstruction for the same acquisition with the addition of
%   transition bands to improve cardiac perfusion.
%--------------------------------------------------------------------------
%   Inputs:      
%       - kSpace_all: Processed kSpace data [nx, nr, nc]
%           - nx: number of measurements along a ray
%           - nr: number of acquired rays
%           - nc: number of PCA coils
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Image: Reconstructed Perfusion Images [nx,ny,nt,nset*nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - Image_PD: Reconstructed Proton Density Images [nx,ny,nt,nset*nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of proton density time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - Image_dia: Reconstructed Diastolic Perfusion Images [nx,ny,nt,nset*nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of diastolic time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - Image_sys: Reconstructed Systolic Perfusion Images [nx,ny,nt,nset*nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of systolic time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Reference:
%       [1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G.
%           Improving ungated steady-state cardiac perfusion using transition
%           bands. Magn Reson Med 2025.
%       [2] Tian Y, Mendes J, Wilson B, Ross A, Ranjan R, DiBella E, Adluru G.
%           Whole-heart, ungated, free-breathing, cardiac-phase-resolved 
%           myocardial perfusion MRI by using Continuous Radial Interleaved
%           simultaneous Multi-slice acquisitions at sPoiled steady-state 
%           (CRIMP). Magn Reson Med 2020;84(6):3071-3087
%--------------------------------------------------------------------------
%   Author:
%       Johnathan Le
%       E-mail: le.johnv@outlook.com
%--------------------------------------------------------------------------

%% Parameter Initialization

[para.Recon.sx,~,para.Recon.no_comp] = size(kSpace_all);

para.Recon.nset = max(para.kSpace_info.set(:))+1;
para.Recon.nSMS = max(para.kSpace_info.phase_mod(:))+1;

nslice = para.Recon.nSMS*para.Recon.nset;
para.Recon.order = 1:para.Recon.nSMS:nslice;
for i=para.Recon.nSMS:-1:2
    para.Recon.order = [para.Recon.order,i:para.Recon.nSMS:nslice];
end
[~,para.Recon.order_back] = sort(para.Recon.order);

%% RING Trajectory Correction using Proton Density Rays
set = para.kSpace_info.set == 0; set = find(set); set = set(1:para.kSpace_info.NumberOfPDlines/para.Recon.nset);
para.trajectory_correction = RING_SMS(kSpace_all(:,set,:), para.kSpace_info.angle_mod(set), para.Recon.nSMS);

%% Proton Density Reconstruction
fprintf('Performing Proton Density Reconstruction...\n')
Image_PD = sliding_window_low_res_recon_for_self_gating_PD(kSpace_all,para);
Image_PD = permute(Image_PD,[1,2,3,5,4]); Image_PD = crop_half_FOV(abs(Image_PD(:,:,:,para.Recon.order)));

%% Cut Proton Density and Non-Steady State Rays
kSpace_all(:,1:para.kSpace_info.NumberOfPDlines+para.Recon.non_steady_state_rays,:) = [];
para.kSpace_info.angle_mod(:,1:para.kSpace_info.NumberOfPDlines+para.Recon.non_steady_state_rays,:) = [];
para.kSpace_info.phase_mod(:,1:para.kSpace_info.NumberOfPDlines+para.Recon.non_steady_state_rays,:) = [];
para.kSpace_info.set(:,1:para.kSpace_info.NumberOfPDlines+para.Recon.non_steady_state_rays,:) = [];
para.kSpace_info.TimeStamp(:,1:para.kSpace_info.NumberOfPDlines+para.Recon.non_steady_state_rays,:) = [];

%% Reference Reconstruction for Self-Gating and Motion Estimation
fprintf('Performing Perfusion Reconstruction...\n')
[Image,para] = sliding_window_low_res_recon_for_self_gating(kSpace_all,para);

%% Self-Gating and Cardiac Phase Fixing
fprintf('Performing Self-Gating and Phase Fixing...\n')
para = get_self_gating_signal(Image,para); [Data_all, para] = fix_Nphase(kSpace_all, Image, para); clear kSpace_all

%% Rigid Translation Estimations using Diastolic Images and Interpolation to All Cardiac Phases
fprintf('Estimating Rigid Translations...\n')
Image_dia = permute(crop_half_FOV(abs(Data_all.first_guess(:,:,para.Recon.bins(para.Recon.Ns2d+1,:),:))),[1,2,4,3]);
[~,shifts_dia] = rigid_reg_3D(Image_dia,para.Recon.self_gating.motion_mask);
shifts = interp_rigid_shifts(shifts_dia,para.Recon.self_gating.resp_signal,para.Recon.bins(para.Recon.Ns2d+1,:));

%% Patch Tracking
fprintf('Performing Patch Tracking...\n')
Data_all = Patch_tracking_3D_with_guide_random_search(Data_all,shifts - mean(shifts,2),para); clearvars -except Data_all Image_PD para

%% Final ADMM Reconstruction
fprintf('Performing ADMM reconstruction...\n')
[Image,para] = LLR_patch_tracking_4D_ADMM(Data_all,para);

%% Save Final Reconstructions
fprintf('Saving Reconstruction Files...\n')
Image = crop_half_FOV(abs(Image)); Image_sys = Image(:,:,para.Recon.bins(1,:),:); Image_dia = Image(:,:,para.Recon.bins(para.Recon.Ns2d+1,:),:);
save(fullfile(para.dir.save_recon_img_dir,para.dir.save_recon_img_name),'Image_PD','Image','Image_sys','Image_dia','para','-v7.3')

return
