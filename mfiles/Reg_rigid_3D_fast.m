function [I2,shifts,step] = Reg_rigid_3D_fast(I1,I0,step,noi,mask,flag)
%--------------------------------------------------------------------------
%   [I2,shifts,step] = Reg_rigid_3D_fast(I1,I0,step,noi,mask,flag)
%--------------------------------------------------------------------------
%   Function to perform 3D image-to-image rigid registration
%--------------------------------------------------------------------------
%   Inputs:      
%       - I1: image to be rigid registered [nx,ny,nset*nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - I0: reference image for registration [nx,ny,nset*nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - step: step size used for each registration update [scalar]
%       - noi: number of iterations used for registration [scalar]
%       - mask: registration mask for region of interest
%       - flag: first iteration flag [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - I2: rigid registered images [nx,ny,nset*nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - shifts: x- and y- rigid shifts used to produce I2 [2,1]
%       - step: step size used after registration update [scalar]
%--------------------------------------------------------------------------

[sx,sy,sz] = size(I0);

a = 30;
s = 3;

kx = cos(2*pi*(0:sx-1)/sx);
ky = cos(2*pi*(0:sy-1)/sy);
kz = cos(2*pi*(0:sz/(sz-1):sz));
W = 2*(kx+ky.'+permute(kz,[3 1 2])-3);
W = (1-a*W).^-s;

N = sum(mask(:));

I2_fft = fftshift3(fft3(fftshift3(I1)));

PhaseX = reshape(-2i*pi*((1:sx)-floor(sx/2)-1)/sx,sx,1);
PhaseY = reshape(-2i*pi*((1:sy)-floor(sy/2)-1)/sy,1,sy);
PhaseZ = reshape(-2i*pi*((1:sz)-floor(sz/2)-1)/sz,1,1,sz);

Phase = exp(PhaseX*0 + PhaseY*0 + PhaseZ*0);
shifts = zeros(3,1);
for iter=1:noi

    I2 = ifftshift3(ifft3(fftshift3(I2_fft.*Phase)));
    I2 = abs(I2);
    
    ddx = 0.5*(I2(3:end,:,:) - I2(1:end-2,:,:));
    ddy = 0.5*(I2(:,3:end,:) - I2(:,1:end-2,:));
    ddz = 0.5*(I2(:,:,3:end) - I2(:,:,1:end-2));
    ddx = cat(1,I2(2,:,:) - I2(1,:,:),ddx,I2(end,:,:) - I2(end-1,:,:));
    ddy = cat(2,I2(:,2,:) - I2(:,1,:),ddy,I2(:,end,:) - I2(:,end-1,:));
    ddz = cat(3,I2(:,:,2) - I2(:,:,1),ddz,I2(:,:,end) - I2(:,:,end-1));
    
    dI = I2-I0;
    
    dx = W.*fft3(ddx.*dI);
    dx = -real(ifft3(dx));

    dy = W.*fft3(ddy.*dI);
    dy = -real(ifft3(dy));
    
    dz = W.*fft3(ddz.*dI);
    dz = -real(ifft3(dz));
    
    dx = dx.*mask;
    dy = dy.*mask;
    dz = dz.*mask;

    dx = -mean(dx(:));
    dy = -mean(dy(:));
    dz = -mean(dz(:))/4;

    d = sos([dx,dy,dz]);
    if d<0.01
        return
    end
    if iter==1 && flag
        step = step/d;
    end
    E(iter) = sum(d(:));
    if iter>1 && E(iter)>E(iter-1)
        step = step/2;
    end
    
    dx = dx*step;
    dy = dy*step;
    dz = dz*step;
    
    
    
    shifts = shifts + [dx;dy;dz];
    Phase = exp(PhaseX*shifts(1) + PhaseY*shifts(2) + PhaseZ*shifts(3));

end
end

function input = fft3(input)

input = fft(input,[],1);
input = fft(input,[],2);
input = fft(input,[],3);
end

function input = ifft3(input)

input = ifft(input,[],1);
input = ifft(input,[],2);
input = ifft(input,[],3);

end

function input = fftshift3(input)

input = fftshift(input,1);
input = fftshift(input,2);
input = fftshift(input,3);
end

function input = ifftshift3(input)

input = ifftshift(input,1);
input = ifftshift(input,2);
input = ifftshift(input,3);
end