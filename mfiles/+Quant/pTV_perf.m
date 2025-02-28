function registered = pTV_perf(Image)
%--------------------------------------------------------------------------
%   registered = pTV_perf(Image)
%--------------------------------------------------------------------------
%   This function performs a group-wise registration of the Perfusion images 
%   using the local nuclear norm metric
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: Unregistered Images [nx,ny,nt]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - registered: Registered Images [nx,ny,nt]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Reference:
%       [1] Vishnevskiy V, Gass T, Szekely G, Tanner C, Goksel O. Isotropic
%           total variation regularization of displacements in parametric image
%           registration. IEEE Trans Med Imaging 2017;36(2):385-395
%--------------------------------------------------------------------------

[nx,ny,nt] = size(Image);

Image = sqrt(abs(double(Image)));

imax = max(max(Image,[],2),[],1);

Image = Image ./ imax;

Image = permute(Image,[1,2,5,4,3]);

%%

opts = [];

opts.pix_resolution = [1.8056,1.8056];
opts.metric = 'local_nuclear';

opts.local_nuclear_patch_size = [4,4];
opts.metric_param = [4,4];
opts.grid_spacing = [4,4]*3;

opts.isoTV = 0.11;
opts.jac_reg = 2;

opts.interp_type = 0;
opts.spline_order = 1;
opts.display = 'off';

opts.border_mask = 20;
opts.k_down = 0.7;

opts.max_iters = 40;

registered = ptv_register(Image,[],opts);

registered = permute(registered,[1,2,5,4,3]);

registered = reshape(registered,[nx,ny,nt]);

registered = gather(single((registered .* imax).^2));

end
