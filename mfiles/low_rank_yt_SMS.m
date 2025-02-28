function update = low_rank_yt_SMS(Image,weight,para)
%--------------------------------------------------------------------------
%   update = low_rank_yt_SMS(Image,weight,para)
%--------------------------------------------------------------------------
%   computes the locally low rank update term
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,1,nsl]
%       - weight: spatial total variation regularization weight [scalar]
%       - para: reconstruction parameters [structure]
%           - Recon
%               - bloc_x: X-direction patch length
%               - bloc_y: Y-direction patch length
%               - tau: value for soft-thresholding
%--------------------------------------------------------------------------
%   Outputs:
%       - update: LLR update term [nx,ny,nt,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------

if weight ~= 0

    [sx,sy,nof,nset,nslice] = size(Image);

    bloc_x = para.Recon.bloc_x;
    bloc_y = para.Recon.bloc_y;
    tau = para.Recon.tau;
    
    update = reshape(Image,[bloc_x,sx/bloc_x,bloc_y,sy/bloc_y,nof,nslice]);
    update = permute(update,[1 3 6 2 4 5]);
    update = reshape(update,[bloc_x*bloc_y*nslice,sx*sy/bloc_x/bloc_y,nof]);
    update = permute(update,[3 1 2]);
    update = gather(update); % CPU is faster here. When using parfor, ~10 times faster than old code.
    
    parfor i=1:sx*sy/bloc_x/bloc_y
        [U,S,V] = svd(update(:,:,i),'econ');
        S = S - tau; S(S < 0) = 0;
        update(:,:,i) = U*S*V';
    end
    
    update = permute(update,[2 3 1]);
    update = reshape(update,[bloc_x,bloc_y,nslice,sx/bloc_x,sy/bloc_y,nof]);
    update = permute(update,[1 4 2 5 6 3]);
    update = reshape(update,[sx,sy,nof,1,nslice]);

    if para.setting.ifGPU
        update = gpuArray(update);
    end
    
else
    update = 0;
end
end
