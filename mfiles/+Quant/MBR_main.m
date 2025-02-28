function [Image_reg,MBI] = MBR_main(Image,numPD,noi)
%--------------------------------------------------------------------------
%   [Image_reg,MBI] = MBR_main(Image,numPD,noi)
%--------------------------------------------------------------------------
%   This function performs an image-to-image registration of the Perfusion
%   images against their corresponding Model-based images by fitting the
%   data to the two-compartment model
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: Unregistered Images [nx,ny,nt]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%       - numPD: number of Proton Density frames in Imaage [scalar]
%       - noi: number of MBR iterations [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - registered: Registered Images [nx,ny,nt]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Reference:
%       [1] Adluru G, Dibella EV, Schabel MC, Model-based registration for
%           dynamic cardiac perfusion MRI. J Magn Reson Imaging
%           2006;24(5):1062-1070
%       [2] Likhite DA, Adluru G, DiBella E, Deformable and Rigid
%           Model-Based Image Registration for Quantitative Cardiac Perfusion.
%           Switzerland: Springer International: 2015.
%--------------------------------------------------------------------------

Image_PD = Image(:,:,1:numPD);
Image_dce = Image(:,:,numPD+1:end);
[sx,sy,nof] = size(Image_dce);

bw = Quant.auto_roi_mbr(Image_dce);

for MBR_iter_no=1:noi

    bldcurve = squeeze(sum(sum(bw.*Image_dce))./sum(sum(bw)));

    tisscurves = reshape(Image_dce,sx*sy,nof);

    nprecontrast = Quant.auto_firstpass(bldcurve'); nprecontrast = 1:nprecontrast;

    init_bld = mean(bldcurve(nprecontrast));

    deltaSI_bldcurve = (bldcurve - init_bld);

    init_tiss = mean(tisscurves(:,nprecontrast),2);
    deltaSI_tisscurve = tisscurves - init_tiss;

    deltaSI_curves = [deltaSI_bldcurve'; deltaSI_tisscurve];

    MBI = Quant.murase_fitg(deltaSI_curves);
    MBI = MBI + init_tiss';
    MBI = reshape(MBI,[nof,sx,sy]);
    MBI = permute(MBI,[2,3,1]);

    Image_dce = Quant.pTV_mbr(Image_dce,MBI);

end

Image_reg = cat(3,Image_PD,Image_dce);

end

