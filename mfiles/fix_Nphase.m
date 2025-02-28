function [Data_all, para] = fix_Nphase(kSpace_all, Image, para)
%--------------------------------------------------------------------------
%   recon_2D_ungated_SPGR_mSMS_ADMM(kSpace_all,para)
%--------------------------------------------------------------------------
%   Function to perform cardiac phase fixing. The number of cardiac phases
%   is fixed for each cardiac cycle based on the extracteed cardiac signal.
%   The preliminary reconstruction and otoher reconstruction parameters are
%   re-grouped to fix the number of frames between each cardiac cycle.
%--------------------------------------------------------------------------
%   Inputs:      
%       - kSpace_all: Processed kSpace data [nx, nr, nc]
%           - nx: number of measurements along a ray
%           - nr: number of acquired rays
%           - nc: number of PCA coils
%       - Image: Preliminary Reconstructed Perfusion Images [nx,ny,nt,nset,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Data_all: Reconstruction variables after cardiac phase fixing [structure]
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
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------

%% pick systole and diastole from cardiac signal
sys = local_max(-para.Recon.self_gating.cardiac_signal);
dia = local_max(para.Recon.self_gating.cardiac_signal);

dia(dia<sys(2)) = []; dia(dia>sys(end-1)) = [];

s2d_correction = dia - sys(2:end-2); s2d_correction = find(s2d_correction < 3);
dia(s2d_correction) = []; sys(s2d_correction+1) = [];

d2s_correction = sys(3:end-1) - dia; d2s_correction = find(d2s_correction < 3);
dia(d2s_correction) = []; sys(d2s_correction+2) = [];

para.Recon.self_gating.sys = sys; para.Recon.self_gating.dia = dia;

para.kSpace_info.TimeStamp = para.kSpace_info.TimeStamp(para.Recon.self_gating.sys(2:end-2));

%% display self-gating result 
figure
plot(para.Recon.self_gating.cardiac_signal, 'LineWidth', 2)
hold on
plot(sys,para.Recon.self_gating.cardiac_signal(sys), '.', 'MarkerSize', 20)
plot(dia,para.Recon.self_gating.cardiac_signal(dia), '.', 'MarkerSize', 20)
legend('cardiac signal', 'self-gated systole', 'self-gated diastole')
xlabel 'Time Frame'
ylabel 'Cardiac Signal [a. u.]'

%% Cardiac Phase Fixing
ray_idx_start_sys = (sys-1)*para.Recon.nor_sl+1;
ray_idx_start_dia = (dia-1)*para.Recon.nor_sl+1;

Ns2d = round(mean((ray_idx_start_dia - ray_idx_start_sys(2:end-2))/para.Recon.nor));
Nd2s = round(mean((ray_idx_start_sys(3:end-1) - ray_idx_start_dia)/para.Recon.nor));
Nphase = Ns2d + Nd2s;

Ncycle = length(ray_idx_start_dia);

tic
for i=1:para.Recon.nset
    set = para.kSpace_info.set == i-1;
    
    kSpace_set = kSpace_all(:,set,:);
    theta_set = para.kSpace_info.angle_mod(set);
    phase_set = para.kSpace_info.phase_mod(set);
    Image_set = squeeze(Image(:,:,:,i,:));
    
    kSpace_reorder = zeros(para.Recon.sx,para.Recon.nor,Nphase,Ncycle,para.Recon.no_comp);
    theta_reorder = zeros(1,para.Recon.nor,Nphase,Ncycle);
    phase_reorder = repmat(0:para.Recon.nSMS-1,[1,para.Recon.nor/para.Recon.nSMS*10,Nphase,Ncycle]);
    Image_interp = zeros(para.Recon.sx,para.Recon.sx,Nphase,Ncycle,para.Recon.nSMS);
    resp_signal_reorder = zeros(Nphase,Ncycle);
    
    for j=1:Ncycle
        sys_start_cycle = ray_idx_start_sys(j+1);
        dia_start_cycle = ray_idx_start_dia(j);
        sys_start_next = ray_idx_start_sys(j+2);
        
        % put the systole rays
        kSpace_reorder(:,1:para.Recon.nor,1,j,:) = kSpace_set(:,sys_start_cycle:sys_start_cycle+para.Recon.nor-1,:);
        theta_reorder(1,1:para.Recon.nor,1,j) = theta_set(sys_start_cycle:sys_start_cycle+para.Recon.nor-1);
        phase_reorder(1,1:para.Recon.nor,1,j) = phase_set(sys_start_cycle:sys_start_cycle+para.Recon.nor-1);
        
        % put the systole to diastole rays
        n_frame_in_between = Ns2d-1;
        nor_in_between = dia_start_cycle - sys_start_cycle - para.Recon.nor + 1;
        nor_in_between_one_frame = ceil(nor_in_between/n_frame_in_between/(para.Recon.nSMS*2))*(para.Recon.nSMS*2);
        
        iphase = 2;
        for k = 1:n_frame_in_between
            ray_start_frame = sys_start_cycle + para.Recon.nor + (k-1)*nor_in_between_one_frame;
            
            kSpace_reorder(:,1:nor_in_between_one_frame,iphase,j,:) = kSpace_set(:,ray_start_frame:ray_start_frame+nor_in_between_one_frame-1,:);
            theta_reorder(1,1:nor_in_between_one_frame,iphase,j) = theta_set(ray_start_frame:ray_start_frame+nor_in_between_one_frame-1);
            phase_reorder(1,1:nor_in_between_one_frame,iphase,j) = phase_set(ray_start_frame:ray_start_frame+nor_in_between_one_frame-1);

            iphase = iphase + 1;
        end
        
        % interp image
        frame = sys(j+1):dia(j);
        Image_for_interp = Image_set(:,:,frame,:);
        Image_for_interp = permute(Image_for_interp,[3,1,2,4]);
        Image_for_interp = Image_for_interp(:,:);
        Image_for_interp = interp1(1:length(frame),Image_for_interp,1:(length(frame)-1)/(Ns2d-1):length(frame));
        Image_for_interp = reshape(Image_for_interp,[Ns2d,para.Recon.sx,para.Recon.sx,para.Recon.nSMS]);
        Image_for_interp = permute(Image_for_interp,[2,3,1,4]);
        Image_interp(:,:,1:Ns2d,j,:) = Image_for_interp;
        
        % interp resp_signal
        resp_signal = para.Recon.self_gating.resp_signal(frame);
        resp_signal_reorder(1:Ns2d,j) = interp1(1:length(frame),resp_signal,1:(length(frame)-1)/(Ns2d-1):length(frame));

        % put the diastole rays
        kSpace_reorder(:,1:para.Recon.nor,iphase,j,:) = kSpace_set(:,dia_start_cycle:dia_start_cycle+para.Recon.nor-1,:);
        theta_reorder(1,1:para.Recon.nor,iphase,j) = theta_set(dia_start_cycle:dia_start_cycle+para.Recon.nor-1);
        phase_reorder(1,1:para.Recon.nor,iphase,j) = phase_set(dia_start_cycle:dia_start_cycle+para.Recon.nor-1);

        iphase = iphase + 1;
        
        % put the diastole to systole rays
        n_frame_in_between = Nd2s-1;
        nor_in_between = sys_start_next - dia_start_cycle - para.Recon.nor + 1;
        nor_in_between_one_frame = ceil(nor_in_between/n_frame_in_between/(para.Recon.nSMS*2))*(para.Recon.nSMS*2);
        
        for k = 1:n_frame_in_between
            ray_start_frame = dia_start_cycle + para.Recon.nor + (k-1)*nor_in_between_one_frame;

            kSpace_reorder(:,1:nor_in_between_one_frame,iphase,j,:) = kSpace_set(:,ray_start_frame:ray_start_frame+nor_in_between_one_frame-1,:);
            theta_reorder(1,1:nor_in_between_one_frame,iphase,j) = theta_set(ray_start_frame:ray_start_frame+nor_in_between_one_frame-1);
            phase_reorder(1,1:nor_in_between_one_frame,iphase,j) = phase_set(ray_start_frame:ray_start_frame+nor_in_between_one_frame-1);
            
            iphase = iphase + 1;
        end
        
        % interp image
        frame = dia(j)+1:sys(j+2)-1;
        Image_for_interp = Image_set(:,:,frame,:);
        Image_for_interp = permute(Image_for_interp,[3,1,2,4]);
        Image_for_interp = Image_for_interp(:,:);
        Image_for_interp = interp1(1:length(frame),Image_for_interp,1:(length(frame)-1)/(Nd2s-1):length(frame));
        Image_for_interp = reshape(Image_for_interp,[Nd2s,para.Recon.sx,para.Recon.sx,para.Recon.nSMS]);
        Image_for_interp = permute(Image_for_interp,[2,3,1,4]);
        Image_interp(:,:,Ns2d+1:Nphase,j,:) = Image_for_interp;
        
        % interp resp_signal
        resp_signal = para.Recon.self_gating.resp_signal(frame);
        resp_signal_reorder(Ns2d+1:Nphase,j) = interp1(1:length(frame),resp_signal,1:(length(frame)-1)/(Nd2s-1):length(frame));

    end
    
    % pre-interpolate
    nor_max = size(kSpace_reorder,2);
    kSpace_reorder = reshape(kSpace_reorder,[para.Recon.sx,nor_max,Nphase*Ncycle,para.Recon.no_comp]);
    theta_reorder = reshape(theta_reorder,[1,nor_max,Nphase*Ncycle]);
    phase_reorder(:,nor_max+1:end,:,:) = [];
    phase_reorder = reshape(phase_reorder,[1,nor_max,Nphase*Ncycle]);
    
    ramp_weight = abs(kSpace_reorder(floor(para.Recon.sx/2)+1,:,:,1)); ramp_weight(ramp_weight ~= 0) = 1;
    Data{i}.ramp_weight = squeeze(sum(ramp_weight,2));
    
    [kx,ky] = get_k_coor(para.Recon.sx,theta_reorder,floor(para.Recon.sx/2)+1);
    
    % RING trajectory correction
    correction = para.trajectory_correction.*permute([cos(theta_reorder);sin(theta_reorder)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);
    
    Data{i}.first_guess = reshape(Image_interp,[para.Recon.sx,para.Recon.sx,Nphase*Ncycle,1,para.Recon.nSMS]);

    Data{i}.kSpace = GROG.SMS_GROG(kSpace_reorder,kx,ky,phase_reorder);
    
    Data{i}.kSpace = orintate_image(Data{i}.kSpace,para.kSpace_info.orintation);

    [Data{i},para] = get_Data_SMS(Data{i},para);
    
end
toc

%% Assignment and Acculumate Slices into a Single Structure

para.Recon.self_gating.resp_signal = resp_signal_reorder(:);
para.Recon.bins = repmat(diag(true(1,Nphase)),[1,Ncycle]);
para.Recon.Ns2d = Ns2d; para.Recon.Nd2s = Nd2s;
para.Recon.Nphase = Nphase; para.Recon.Ncycle = Ncycle;

Data_all = Data{1};
for i=2:para.Recon.nset
    Data_all.kSpace = cat(6,Data_all.kSpace,Data{i}.kSpace);
    Data_all.mask = cat(6,Data_all.mask,Data{i}.mask);
    Data_all.sens_map = cat(6,Data_all.sens_map,Data{i}.sens_map);
    Data_all.first_est = cat(6,Data_all.first_est,Data{i}.first_est);
    Data_all.first_guess = cat(6,Data_all.first_guess,Data{i}.first_guess);
    if isfield(Data_all,'ramp_filter')
        Data_all.ramp_filter = cat(6,Data_all.ramp_filter,Data{i}.ramp_filter);
    end
    if isfield(Data_all,'sms_filter')
        Data_all.sms_filter = cat(6,Data_all.sms_filter,Data{i}.sms_filter);
    end
end
Data_all.first_guess = single(Data_all.first_guess(:,:,:,para.Recon.order)); para.Recon.scale = max(abs(vec(Data_all.first_est)));

Data_all = rmfield(Data_all,'first_est'); Data_all = rmfield(Data_all,'ramp_weight');


