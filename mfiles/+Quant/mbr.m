function [Image_reg,MBI] = mbr(Image,numPD)
%--------------------------------------------------------------------------
%   [Image_reg,MBI] = mbr(Image,numPD)
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

if ~exist('numPD')
    numPD = 1;
end
noi=2;  % total number of MBR iterations

[Image_reg,MBI] = Quant.MBR_main(Image,numPD,noi);
