function [Image, para] = sliding_window_low_res_recon_for_self_gating(kSpace_all,para)
%--------------------------------------------------------------------------
%   [Image, para] = sliding_window_low_res_recon_for_self_gating(kSpace_all,para)
%--------------------------------------------------------------------------
%   Function to perform preliminary STCR reconstruction prior to the final
%   patch-based total variation and LLR regularization ADMM reconstruction.
%   Used for self-gating, respiratory binning, and motion compensation
%   prior to final reconstruction.
%--------------------------------------------------------------------------
%   Inputs:      
%       - kSpace_all: Processed kSpace data [nx, nr, nc]
%           - nx: number of measurements along a ray
%           - nr: number of acquired rays
%           - nc: number of PCA coils
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Image: Preliminary Reconstructed Perfusion Images [nx,ny,nt,nset,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------

%% Pre-calculate GROG gridding operators for interpolating radial data to Cartesian
for i=1:para.Recon.nset

    set = para.kSpace_info.set == i-1;
    
    kSpace_set = kSpace_all(:,set,:);
    theta_set = para.kSpace_info.angle_mod(set);
    phase_set = para.kSpace_info.phase_mod(set);
    
    for j=1:para.Recon.nSMS
        phase_SMS = phase_set == j-1;
        
        kSpace_SMS = kSpace_set(:,phase_SMS,:);
        kSpace_SMS = permute(kSpace_SMS,[1,2,4,3]);
        theta_SMS = theta_set(phase_SMS);

        [kx, ky] = get_k_coor(para.Recon.sx,theta_SMS,floor(para.Recon.sx/2)+1);

        if isfield(para,'trajectory_correction')
            correction = para.trajectory_correction.*permute([cos(theta_SMS);sin(theta_SMS)],[4,1,2,3]);
            correction = squeeze(sum(correction,2));
            kx = kx - correction(1,:,:);
            ky = ky - correction(2,:,:);
        end

        [para.G{i,j}.Gx, para.G{i,j}.Gy] = GROG.get_Gx_Gy(kSpace_SMS, kx, ky);
    end
end

%% Data binning and radial to Cartesian interpolation
for i=1:para.Recon.nset
    set = para.kSpace_info.set == i-1;

    kSpace_set = kSpace_all(:,set,:);
    theta_set = para.kSpace_info.angle_mod(set);
    phase_set = para.kSpace_info.phase_mod(set);
    
    nor_total = size(kSpace_set,2);
    
    nof = floor((nor_total-(para.Recon.nor-para.Recon.nor_sl))/para.Recon.nor_sl);
    nor_total = nof*para.Recon.nor_sl+(para.Recon.nor-para.Recon.nor_sl);
    
    kSpace_set(:,nor_total+1:end,:) = [];
    theta_set(nor_total+1:end) = [];
    phase_set(nor_total+1:end) = [];
    
    kSpace = zeros(para.Recon.sx,para.Recon.nor,nof,para.Recon.no_comp);
    theta = zeros(1,para.Recon.nor,nof);
    phase = zeros(1,para.Recon.nor,nof);

    for j=1:nof
        ray = (j-1)*para.Recon.nor_sl+1:(j-1)*para.Recon.nor_sl+para.Recon.nor;
        kSpace(:,:,j,:) = kSpace_set(:,ray,:);
        theta(1,:,j) = theta_set(ray);
        phase(1,:,j) = phase_set(ray);
    end
    
    ramp_weight = abs(kSpace(floor(para.Recon.sx/2)+1,:,:,1)); ramp_weight(ramp_weight ~= 0) = 1;
    Data{i}.ramp_weight = squeeze(sum(ramp_weight,2));
    
    [kx,ky] = get_k_coor(para.Recon.sx,theta,floor(para.Recon.sx/2)+1);
    
    if isfield(para,'trajectory_correction')
        correction = para.trajectory_correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
        correction = squeeze(sum(correction,2));
        kx = kx - correction(1,:,:);
        ky = ky - correction(2,:,:);
    end
    
    Data{i}.kSpace = GROG.SMS_GROG(kSpace,kx,ky,phase,para.G(i,:));
    
    Data{i}.kSpace = orintate_image(Data{i}.kSpace,para.kSpace_info.orintation);
    
    if para.Recon.nSMS == 1
        [Data{i},para] = get_Data(Data{i},para);
    else
        [Data{i},para] = get_Data_SMS(Data{i},para);
    end

end
clearvars -except Data para

para.kSpace_info.TimeStamp = para.kSpace_info.TimeStamp(1:para.Recon.nset*para.Recon.nor_sl:end);

%% Preliminary STCR reconstruction

% Preliminary STCR reconstruction replaces the initial zero-filled image
% estimate prior to the final patch-based total variation and LLR 
% regularization ADMM reconstruction.
siz = size(Data{1}.first_est); siz(4) = para.Recon.nset; Image = zeros(siz,'single');
for i=1:para.Recon.nset
    scale = max(abs(Data{i}.first_est(:)));
    para.Recon.weight_tTV = scale*para.weight_tTV;
    para.Recon.weight_sTV = scale*para.weight_sTV;
    [Image(:,:,:,i,:),para] = STCR_conjugate_gradient(Data{i},para);
end
