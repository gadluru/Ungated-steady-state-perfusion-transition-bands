function registered = pTV_PD_reg(Image)
%--------------------------------------------------------------------------
%   registered = pTV_PD_reg(Image)
%--------------------------------------------------------------------------
%   This function performs a group-wise registration of the Proton Density 
%   images using the local nuclear norm metric
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

Image = sqrt(abs(double(Image)));

Image = permute(Image,[1,2,5,4,3]);

imax = max(max(Image,[],2),[],1);

Image = Image ./ imax;

%%

opts = [];

opts.pix_resolution = [1.8056,1.8056];
opts.metric = 'local_nuclear';

opts.local_nuclear_patch_size = [4,4];
opts.metric_param = [4,4];
opts.grid_spacing = [4,4]*3;

opts.isoTV = 3e-2;
opts.jac_reg = 2;

opts.interp_type = 0;
opts.spline_order = 1;
opts.display = 'off';

opts.border_mask = 20;
opts.k_down = 0.7;

opts.max_iters = 40;

%%
registered = ptv_register(Image,[],opts);

registered = gather(single(mean((registered .* imax).^2,5)));

end
