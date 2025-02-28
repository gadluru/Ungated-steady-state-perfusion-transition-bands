function [kSpace_cart,G] = SMS_GROG(kSpace_all, kx, ky, phase_mod, G_pre_calculated)
%--------------------------------------------------------------------------
%   kSpace_cart = SMS_GROG(kSpace_all, kx, ky, phase_mod, G_pre_calculated)
%--------------------------------------------------------------------------
%   GROG interpolation function for converison of 2D MRI radial k-space
%   data to 2D MRI cartesian k-space data
%--------------------------------------------------------------------------
%   Inputs:      
%       - kSpace_all: Processed kSpace data [nx,nr,nt,nc]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame
%           - nt: number of frames
%           - nc: number of PCA coils
%       - kx: k-space coordinates kx, normalized to 2x matrix size [nx,nr,nt]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame
%           - nt: number of frames
%       - ky: k-space coordinates ky, normalized to 2x matrix size [nx,nr,nt]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame
%           - nt: number of frames
%       - phase_mod: k-space coordinates ky, normalized to 2x matrix size [1,nr,nt]
%           - nr: number of rays per frame
%           - nt: number of frames
%       - G_pre_calculated: GROG gridding operators [structure]
%--------------------------------------------------------------------------
%   Output:
%       - kSpace_cart: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nc: number of PCA coils
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------
%   Reference:
%       - Ye Tian, et al. (2019) PLOS ONE 14(2): e0211738.
%       - Nicole Seiberlich, et al. (2008) MRM 59:930-935.
%--------------------------------------------------------------------------
% Copyright: University of Utah Cardiovascular MRI Group
% https://medicine.utah.edu/radiology/radiology-research/research-labs/dibella2/
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

nSMS = max(phase_mod(:))+1;

[sx,nor,nof,nc,~,NSlice] = size(kSpace_all);

G = cell(1,nSMS);

for i=1:NSlice
    
    kSpace_slice = kSpace_all(:,:,:,:,i);
    
    for j=1:nSMS
        SMS = phase_mod == j-1;
        
        kSpace_SMS = kSpace_slice(repmat(SMS,[sx 1 1 nc]));
        kSpace_SMS = reshape(kSpace_SMS,[sx,nor/nSMS,nof,nc]);
        
        kx_SMS = kx(repmat(SMS,[sx 1 1]));
        ky_SMS = ky(repmat(SMS,[sx 1 1]));
        
        kx_SMS = reshape(kx_SMS,[sx nor/nSMS nof]);
        ky_SMS = reshape(ky_SMS,[sx nor/nSMS nof]);
        
        mask = abs(kSpace_SMS(floor(sx/2)+1,:,:,1)); mask(mask ~= 0) = 1; mask = squeeze(sum(mask,2))';
        
        for k=unique(mask)
            idx1 = (mask == k);
            
            kSpace_block = kSpace_SMS(:,:,idx1,:);
            kx_block = kx_SMS(:,:,idx1);
            ky_block = ky_SMS(:,:,idx1);

            idx2 = squeeze(sum(sum(sum(abs(kSpace_block),1),3),4)) ~= 0;

            kSpace_block = kSpace_block(:,idx2,:,:);
            kx_block = kx_block(:,idx2,:);
            ky_block = ky_block(:,idx2,:);

            if exist('G_pre_calculated')
                G{j} = G_pre_calculated{j};
                kSpace_cart(:,:,idx1,:,:,:,j) = GROG.GROG_Dictionary_interp_new(kSpace_block,G_pre_calculated{j}.Gx,G_pre_calculated{j}.Gy,kx_block,ky_block);
            else
                [G{j}.Gx, G{j}.Gy] = GROG.get_Gx_Gy(kSpace_block, kx_block, ky_block);
                kSpace_cart(:,:,idx1,:,:,:,j) = GROG.GROG_Dictionary_interp_new(kSpace_block,G{j}.Gx,G{j}.Gy,kx_block,ky_block);
            end
        end
    end
end

sx_over = size(kSpace_cart,1);
kSpace_cart = reshape(kSpace_cart,[sx_over,sx_over,nof,nc,1,NSlice,nSMS]);
if size(kSpace_cart,1)~=sx
    kSpace_cart([1,end],:,:,:,:,:,:) = [];
    kSpace_cart(:,[1,end],:,:,:,:,:) = [];
end
