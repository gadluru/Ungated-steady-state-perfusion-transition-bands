function Data_all = Patch_tracking_3D_with_guide_random_search(Data_all,offset_init,para)
%--------------------------------------------------------------------------
%   Data_all = Patch_tracking_3D_with_guide_random_search(Data_all,offset_init,para)
%--------------------------------------------------------------------------
%   This function performs 3D patch tracking for each cardiac phase along
%   cardiac cycles, where patches with significant motion or dynamic
%   changes are tracked (based on temporal standard deviation maps). The
%   interpolated diastolic translations are used to reduce the search
%   window of each patch.
%--------------------------------------------------------------------------
%   Inputs:      
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
%       - offset_init: interpolated diastolic rigid translations [3,nt]
%           - nt: number of time frames after cardiac phase fixing
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
%--------------------------------------------------------------------------

tic

patch_size = [5,5,2]; search_size = [9,9,4]; patch_shift = [5,5,2];

[nx,ny,nt,nsl] = size(Data_all.first_guess);

Nphase = size(para.Recon.bins,1); Ncycle = nt/Nphase;

im_bin = reshape(Data_all.first_guess,[nx,ny,Nphase,Ncycle,nsl]);
im_bin = permute(im_bin,[1,2,5,3,4]);

x_begin_all = 1:patch_shift(1):nx-patch_size(1)+1;
y_begin_all = 1:patch_shift(2):ny-patch_size(2)+1;
z_begin_all = 1:patch_shift(3):nsl-patch_size(3)+1;

if z_begin_all(end) + patch_size(3)-1 < nsl
    z_begin_all(end+1) = nsl-patch_size(3) + 1;
end
if x_begin_all(end) + patch_size(1)-1 < nx
    x_begin_all(end+1) = nx-patch_size(1) + 1;
end
if y_begin_all(end) + patch_size(2)-1 < ny
    y_begin_all(end+1) = ny-patch_size(2) + 1;
end
Nx = length(x_begin_all); Ny = length(y_begin_all); Nz = length(z_begin_all);

std_map = std(im_bin,1,5); int_map = sum(abs(im_bin),5);

offset_init = reshape(offset_init,[3,Nphase,Ncycle]);
offset_init = offset_init - offset_init(:,:,1);
offset_init(:,:,2:end) = diff(offset_init,1,3);
offset_init = round(offset_init);
offset_init = -offset_init;

%% Extract possible motion-tracked patches and their linear indices
idx_all = zeros([patch_size,Nx,Ny,Nz,Nphase]);
int_patch = zeros([patch_size,Nx,Ny,Nz,Nphase]);
std_patch = zeros([patch_size,Nx,Ny,Nz,Nphase]);
parfor iphase = 1:Nphase
    std_phase = std_map(:,:,:,iphase);
    int_phase = int_map(:,:,:,iphase);
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                x = x_begin_all(i):x_begin_all(i)+patch_size(1)-1;
                y = y_begin_all(j):y_begin_all(j)+patch_size(2)-1;
                z = z_begin_all(k):z_begin_all(k)+patch_size(3)-1;
                
                [x,y,z] = ndgrid(x,y,z);
                idx = sub2ind([nx,ny,nsl],x,y,z);
                
                idx_all(:,:,:,i,j,k,iphase) = idx;
                int_patch(:,:,:,i,j,k,iphase) = int_phase(idx);
                std_patch(:,:,:,i,j,k,iphase) = std_phase(idx);
            end
        end
    end
end
idx_all = reshape(idx_all,[patch_size,Nx*Ny*Nz,Nphase]);
int_patch = reshape(int_patch,[patch_size,Nx*Ny*Nz,Nphase]);
std_patch = reshape(std_patch,[patch_size,Nx*Ny*Nz,Nphase]);

%% Patch Tracking Algorithm
offset_all = {}; phase_all = zeros(1,Nphase); mask_all = zeros(nx,ny,nsl,Nphase,Ncycle);
parfor iphase = 1:Nphase

    Image_phase = squeeze(im_bin(:,:,:,iphase,:));
    offset_phase = squeeze(offset_init(:,iphase,:))';
    idx_all_phase = idx_all(:,:,:,:,iphase);
    int_patch_phase = int_patch(:,:,:,:,iphase);
    std_patch_phase = std_patch(:,:,:,:,iphase);
    
    idx = sum(sum(sum(abs(int_patch_phase))));
    idx_std = sum(sum(sum(abs(std_patch_phase))));
    keep = (idx_std > max(idx_std)/5) & (idx > max(idx)/5);
    
    idx_all_phase = idx_all_phase(:,:,:,keep);
    idx_add = (0:Ncycle-1)*nx*ny*nsl;
    idx_all_phase = idx_all_phase + permute(idx_add,[1,3,4,5,2]);
    
    Npatch = sum(keep(:));
    idx_all_phase = reshape(idx_all_phase,[patch_size,Npatch*Ncycle]);
    
    [~,searching_order] = sort(rand([1,Npatch*Ncycle]));
    
    %% patch search
    
    n = 0; mask = zeros(size(Image_phase)); offset_all_phase = zeros(3,Ncycle,Npatch*Ncycle);
    for i=1:Npatch*Ncycle

        order = searching_order(i);
        idx_patch = idx_all_phase(:,:,:,order);
        flag = sum(vec(mask(idx_patch) == 0));
        
        if flag
            n = n+1;

            [x,y,z,c] = ind2sub([nx,ny,nsl,Ncycle],idx_patch);

            x_begin = x(1); y_begin = y(1); z_begin = z(1); c = c(1);
            offset_all_phase(:,c,n) = [x_begin,y_begin,z_begin];

            mask(idx_patch) = mask(idx_patch) + 1;
            if c~=1
                for icycle=c-1:-1:1
                    offset_cycle = -offset_phase(icycle+1,:);
                    
                    [x_window_begin,x_window_end,y_window_begin,y_window_end,z_window_begin,z_window_end] = get_search_window(x_begin,y_begin,z_begin,search_size,patch_size,offset_cycle,nx,ny,nsl);
                    
                    patch = Image_phase(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle+1);
                    search_window = Image_phase(x_window_begin:x_window_end,y_window_begin:y_window_end,z_window_begin:z_window_end,icycle);
                    
                    offset = patch_search_one_3D(patch,search_window);
                    
                    x_begin = x_window_begin-1+offset(1); y_begin = y_window_begin-1+offset(2); z_begin = z_window_begin-1+offset(3);
                    offset_all_phase(:,icycle,n) = [x_begin,y_begin,z_begin];

                    mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle) = mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle) + 1;
                end
            end
            if c~=Ncycle
                x_begin = x(1); y_begin = y(1); z_begin = z(1);
                for icycle=c:1:Ncycle-1
                    offset_cycle = offset_phase(icycle+1,:);
                    
                    [x_window_begin,x_window_end,y_window_begin,y_window_end,z_window_begin,z_window_end] = get_search_window(x_begin,y_begin,z_begin,search_size,patch_size,offset_cycle,nx,ny,nsl);
                    
                    patch = Image_phase(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle);
                    search_window = Image_phase(x_window_begin:x_window_end,y_window_begin:y_window_end,z_window_begin:z_window_end,icycle+1);
                    
                    offset = patch_search_one_3D(patch,search_window);
                    
                    x_begin = x_window_begin-1+offset(1); y_begin = y_window_begin-1+offset(2); z_begin = z_window_begin-1+offset(3);
                    offset_all_phase(:,icycle+1,n) = [x_begin,y_begin,z_begin];

                    mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle+1) = mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle+1) + 1;
                end
            end
        end
    end
    
    idx = sum(reshape(offset_all_phase,[3*Ncycle,Npatch*Ncycle]),1); idx = find(idx == 0); idx = idx(1) - 1;
    offset_all_phase = offset_all_phase(:,:,1:idx);
    
    phase_all(iphase) = n; offset_all{iphase} = offset_all_phase; mask_all(:,:,:,iphase,:) = mask;
end
offset_all = cat(3,offset_all{:}); Npatch = sum(phase_all);

%% Get linear indices of motion-tracked patches

phase_idx = [];
for i=1:Nphase
    phase_idx = cat(2,phase_idx,ones(1,phase_all(i))*i);
end

patch_begin_all = offset_all; patch_end_all = offset_all + reshape(patch_size(1:3),3,1) - 1;

idx = zeros([patch_size,Ncycle,Npatch]);
parfor i=1:Npatch
    for j=1:Ncycle
        x = patch_begin_all(1,j,i):patch_end_all(1,j,i);
        y = patch_begin_all(2,j,i):patch_end_all(2,j,i);
        z = patch_begin_all(3,j,i):patch_end_all(3,j,i);
        p = phase_idx(i);

        [x,y,z,p] = ndgrid(x,y,z,p);

        idx(:,:,:,j,i) = sub2ind([nx,ny,nsl,Nphase],x,y,z,p);
    end
end
idx = reshape(idx,prod(patch_size),Ncycle,Npatch);
idx_add = nx*ny*nsl*Nphase*(0:Ncycle-1); idx = idx+idx_add;

mask = reshape(mask_all,[nx,ny,nsl,Nphase*Ncycle]); mask_LR = mask;
mask_LR(mask_LR==0) = 1; mask_LR = 1./mask_LR;

Data_all.llr.MT.Npatch = Npatch;
Data_all.llr.MT.idx = int32(idx);
Data_all.llr.MT.mask_intensity = single(mask_LR);
Data_all.llr.MT.mask = logical(mask);

%% Get linear indices of non-moving patches

phase_idx = [];
for i=1:Nphase
    phase_idx = cat(2,phase_idx,ones(1,Nx*Ny*Nz)*i);
end

Npatch = length(phase_idx);

patch_begin_all = zeros(3,1,Nx,Ny,Nz);
parfor i=1:Nx
    for j=1:Ny
        for k=1:Nz
            patch_begin_all(:,1,i,j,k) = cat(1,x_begin_all(i),y_begin_all(j),z_begin_all(k));
        end
    end
end
patch_begin_all = repmat(reshape(repmat(patch_begin_all,[1,Ncycle,1,1,1]),[3,Ncycle,Nx*Ny*Nz]),[1,1,Nphase]);
patch_end_all = patch_begin_all + reshape(patch_size(1:3),3,1) - 1;

idx = zeros([patch_size,Ncycle,Npatch]);
parfor i=1:Npatch
    for j=1:Ncycle
        x = patch_begin_all(1,j,i):patch_end_all(1,j,i);
        y = patch_begin_all(2,j,i):patch_end_all(2,j,i);
        z = patch_begin_all(3,j,i):patch_end_all(3,j,i);
        p = phase_idx(i);

        [x,y,z,p] = ndgrid(x,y,z,p);

        idx(:,:,:,j,i) = sub2ind([nx,ny,nsl,Nphase],x,y,z,p);
    end
end
idx = reshape(idx,prod(patch_size),Ncycle,Npatch);
idx_add = nx*ny*nsl*Nphase*(0:Ncycle-1); idx = idx+idx_add;

mask = accumarray(vec(idx),ones(size(vec(idx))),[numel(Data_all.first_guess),1]);
mask = reshape(mask,[nx,ny,nsl,Nphase*Ncycle]);
mask_LR = mask; mask_LR(mask_LR==0) = 1; mask_LR = 1./mask_LR;

Data_all.llr.NM.Npatch = Npatch;
Data_all.llr.NM.idx = int32(idx);
Data_all.llr.NM.mask_intensity = single(mask_LR);
Data_all.llr.NM.mask = logical(mask);

Data_all.size.Nphase = Nphase; Data_all.size.Ncycle = Ncycle;

Data_all.size.order = para.Recon.order; Data_all.size.order_back = para.Recon.order_back;

Data_all.size.nSMS = para.Recon.nSMS; Data_all.size.nset = para.Recon.nset;

toc

end










