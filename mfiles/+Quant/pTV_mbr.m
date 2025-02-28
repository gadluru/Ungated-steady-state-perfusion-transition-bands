function registered = pTV_mbr(Image,MBI)
%--------------------------------------------------------------------------
%   registered = pTV_mbr(Image,MBI)
%--------------------------------------------------------------------------
%   This function performs image-to-image registration of the Perfusion images 
%   using their corresponding model-based images as a reference.
%   Registration is performed using a local cross correlation metric.
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: Unregistered Images [nx,ny,nt]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%       - MBI: Model-based reference images [nx,ny,nt]
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

Image = sqrt(abs(double(Image))); MBI = sqrt(abs(double(MBI)));

imax = max(max(Image,[],2),[],1); mmax = max(max(MBI,[],2),[],1);

Image = Image ./ imax; MBI = MBI ./ mmax;

%%

opts = [];

opts.pix_resolution = [1.8056,1.8056];
opts.metric = 'loc_cc_fftn';

opts.loc_cc_approximate = true;
opts.loc_cc_abs = true;

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

%%

registered = zeros(nx,ny,nt);
for i=1:nt
    registered(:,:,i) = ptv_register(Image(:,:,i),MBI(:,:,i),opts);
end

registered = gather(single((registered.*imax).^2));

end

