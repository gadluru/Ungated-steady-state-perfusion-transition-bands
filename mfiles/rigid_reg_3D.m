function [Image_reg,shifts] = rigid_reg_3D(Image,mask)
%--------------------------------------------------------------------------
%   [registered,shifts] = rigid_reg_IR_GS(Image,MBI,para)
%--------------------------------------------------------------------------
%   Function to perform 3D rigid registration of the image series
%   interpolated along the slice dimension to improve registration. Each 
%   time frame of the 3D volume is registered to its following 3D volume.
%   The reference 3D volume is updated by averaging with the previously
%   registered 3D volumes.
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: image series to be rigid registered [nx,ny,nt,nset*nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - mask: registration mask for cost update [nx,ny,nset*nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%--------------------------------------------------------------------------
%   Outputs:
%       - Image_reg: rigid registered images [nx,ny,nt,nset*nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%       - shifts: calculated x- and y- rigid shifts [3,noi,nt]
%           - noi: number of iterations
%           - nt: number of time frames
%--------------------------------------------------------------------------

[nx,ny,nsl,nt] = size(Image);

Image = cat(3,Image(:,:,[1,1],:),Image,Image(:,:,[end,end,end],:));
Image = interp_2_slices_in_between(Image); Image_reg = Image;

mask = cat(3,zeros(nx,ny,2),mask,zeros(nx,ny,3));
mask = interp_2_slices_in_between(mask);
mask(:,:,[1:2,end-2:end]) = false;

noi = 4; step = 2; first_iter_flag = 1;
shifts = zeros(3,noi,nt); step = ones(1,nt)*step;

tic
for iter_no = 1:noi
    Image_ref = Image_reg(:,:,:,end);
    for i=nt-1:-1:1
        tic;
        [Image_reg(:,:,:,i),shifts(:,iter_no,i),step(i)] = Reg_rigid_3D_fast(Image(:,:,:,i),Image_ref,step(i),10,mask,first_iter_flag);
        Image_ref = (Image_ref*(nt-i) + Image_reg(:,:,:,i))/(nt-i+1);
    end
    Image = Image_reg;
    first_iter_flag = 0;
end
toc

Image_reg = Image_reg(:,:,1:3:end,:);
Image_reg = Image_reg(:,:,3:end-3,:);

end

function Image = interp_2_slices_in_between(Image)
[sx,sy,sz,nof] = size(Image);
if nof>1
    [x,y,z,f] = ndgrid(1:sx,1:sy,1:sz,1:nof);
    [xv,yv,zv,fv] = ndgrid(1:sx,1:sy,1:1/3:sz,1:nof);
    Image = interpn(x,y,z,f,Image,xv,yv,zv,fv);
else
    [x,y,z] = ndgrid(1:sx,1:sy,1:sz);
    [xv,yv,zv] = ndgrid(1:sx,1:sy,1:1/3:sz);
    Image = interpn(x,y,z,Image,xv,yv,zv);
end

end