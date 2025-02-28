
function para = get_self_gating_signal(Image,para)
%--------------------------------------------------------------------------
%   para = get_self_gating_signal(Image,para)
%--------------------------------------------------------------------------
%   Function to perform self-gating and to obtain the respiration signal
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: Preliminary Reconstructed Perfusion Images [nx,ny,nt,nset,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - para: reconstruction parameters [structure]
%           - Recon.self_gating
%               - cardiac_signal: extracted cardiac signal [nt,1]
%               - resp_signal: extracted respiration signal [nt,1]
%               - cardiac_mask: extracted cardiac mask [nx,ny,nset*nsl]
%               - motion_mask: extracted motion mask [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%--------------------------------------------------------------------------

%% Generate Motion Mask and Cardiac Mask

Image = permute(Image,[1,2,3,5,4]);
Image = Image(:,:,:,para.Recon.order);

image_for_std = abs(crop_half_FOV(Image(:,:,1:28,:)));

[nx,ny,nt,nsl] = size(image_for_std);

image_for_std = reshape(image_for_std,[nx,ny,7,4,nsl]);

std_map = squeeze(sum(std(image_for_std,0,3),4));
std_map = imgaussfilt(std_map,5); std_map = std_map.*fspecial('gaussian',nx,nx/10);

mask = zeros(nx,ny,nsl);
for i=1:nsl
    [x,y] = find(std_map(:,:,i)==max(max(std_map(:,:,i))));
    mask(x,y,i) = 1; mask(:,:,i) = bwdist(mask(:,:,i));
end
motion_mask = mask<35; cardiac_mask = mask<25;

%% Get Cardiac Signal

cardim = abs(crop_half_FOV(Image));

[nx,ny,nt,nsl] = size(cardim);

dt = para.Recon.nor_sl * (para.kSpace_info.Protocol.alTR(2)/1000) * para.Recon.nset;
Fs = 1000/dt;
df = Fs/nt;

cardiac_signal = cardim .* permute(cardiac_mask,[1,2,4,3]);
cardiac_signal = permute(cardiac_signal,[1,2,4,3]);
cardiac_signal = reshape(cardiac_signal,[nx*ny*nsl,nt]);
cardiac_signal = bandpass(sum(cardiac_signal,1),[0.5,2.2],Fs)';

%% Get Respiration Signal

[nx,ny,nt,nsl] = size(Image);

resp_range = round((0.2/df):(0.5/df));

if resp_range==0
    resp_signal = ones(1,nt);
else
    resp_signal = squeeze(sum(abs(Image),2));
    resp_signal = permute(resp_signal,[1,3,2]);
    resp_signal = reshape(resp_signal,[nx*nsl,nt]);
    resp_coeff = pca(resp_signal);
    resp_signal_fft = fft(resp_coeff(:,1:10),[],1);

    [~,idx] = max(max(abs(resp_signal_fft(resp_range,:))),[],2);

    [~,idx_max] = max(abs(resp_signal_fft(resp_range,idx)),[],1);
    peak_frequency = (resp_range(1) + idx_max - 1)*df/Fs;

    resp_signal = lowpass(resp_coeff(:,idx),peak_frequency);
end

%% Assignment to Structure

para.Recon.self_gating.cardiac_signal = imgaussfilt(cardiac_signal,2);
para.Recon.self_gating.resp_signal = imgaussfilt(resp_signal,2);

para.Recon.self_gating.motion_mask = motion_mask;
para.Recon.self_gating.cardiac_mask = cardiac_mask;

end