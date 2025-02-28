function ramp_filter = ramp_filter_for_pre_interp(Data,para)
%--------------------------------------------------------------------------
%   filter = ramp_IR_filter_for_pre_interp(para)
%--------------------------------------------------------------------------
%   Generates ramp filter to compensate for oversampling of k-space center
%   due to radial acquisition
%--------------------------------------------------------------------------
%   Inputs:      
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - filter: radial sampling filter [nx,ny,nt]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%--------------------------------------------------------------------------

[X,Y] = meshgrid(1:para.Recon.kSpace_size(1),1:para.Recon.kSpace_size(2));
X = X - (para.Recon.kSpace_size(1)+1)/2;
Y = Y - (para.Recon.kSpace_size(2)+1)/2;
ramp_filter = repmat(sqrt(X.^2 + Y.^2),[1,1,length(Data.ramp_weight)]);
fully_sampled_radius = (Data.ramp_weight)*2/pi;
for i=1:length(Data.ramp_weight)
    filter = ramp_filter(:,:,i);
    filter(filter < fully_sampled_radius(i)) = fully_sampled_radius(i);
    ramp_filter(:,:,i) = filter;
end
ramp_filter = ramp_filter/(para.Recon.sx+1)*2;
ramp_filter(ramp_filter>1) = 1;
ramp_filter = fftshift2(ramp_filter);
ramp_filter = single(ramp_filter);

end