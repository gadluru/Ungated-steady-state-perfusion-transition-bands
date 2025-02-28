function Q = get_mask_AHA(Q)
%--------------------------------------------------------------------------
%   Q = get_mask_AHA(Q)
%--------------------------------------------------------------------------
%   This function performs manual AHA segmentation of each slice 
%--------------------------------------------------------------------------
%   Inputs:      
%       - Q: Quantification parameters [struct]
%           -Tissue
%               -Image_PD: Proton Density images [nx,ny,nset*nsl]
%               -Image_MBR: Perfusion images [nx,ny,nset*nsl,nt]
%               -SI_curve_second_all: raw signal Perfusion images [nx,ny,nset*nsl,nt]
%               -T1_curve_second_all: T1 images corresponding to Image_MBR [nx,ny,nset*nsl,nt]
%               -Gd_curve_second_all: [Gd] images corresponding to Image_MBR [nx,ny,nset*nsl,nt]
%               -Mask_heart: Mask for blood flow quantification [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%           -TCM
%               -para_low: lower bound of 2CM model parameters [1,4]
%               -para_up: upper bound of 2CM model parameters [1,4]
%               -para_init: initial estimate of 2CM model parameters [1,4]
%               -bld_flow_map: blood flow parameter map (identical to ktrans) [nx,ny,nset*nsl]
%               -kep_map: kep parameter map [nx,ny,nset*nsl]
%               -vb_map: vb parameter map [nx,ny,nset*nsl]
%               -dt_map: dt parameter map [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%--------------------------------------------------------------------------
%   Outputs:
%       - Q: Quantification parameters [struct]
%           -Tissue
%               -Image_PD: Proton Density images [nx,ny,nset*nsl]
%               -Image_MBR: Perfusion images [nx,ny,nset*nsl,nt]
%               -SI_curve_second_all: raw signal Perfusion images [nx,ny,nset*nsl,nt]
%               -T1_curve_second_all: T1 images corresponding to Image_MBR [nx,ny,nset*nsl,nt]
%               -Gd_curve_second_all: [Gd] images corresponding to Image_MBR [nx,ny,nset*nsl,nt]
%               -Mask_heart: Mask for blood flow quantification [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%           -TCM
%               -para_low: lower bound of 2CM model parameters [1,4]
%               -para_up: upper bound of 2CM model parameters [1,4]
%               -para_init: initial estimate of 2CM model parameters [1,4]
%               -bld_flow_map: blood flow parameter map (identical to ktrans) [nx,ny,nset*nsl]
%               -kep_map: kep parameter map [nx,ny,nset*nsl]
%               -vb_map: vb parameter map [nx,ny,nset*nsl]
%               -dt_map: dt parameter map [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%           -mask
%               -Seg_mask: AHA regional binary mask for each slice [cell]
%               -Mask_all: myocardial binary mask for each slice [nx,ny,nset*nsl]
%               -Image_show: Display image for Bullseye [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%--------------------------------------------------------------------------

bld_flow_map = Q.TCM.bld_flow_map;
Nslice = size(Q.Tissue.Image_MBR,3);

fprintf('Draw LV: \n')
curve = Quant.compare_curve_same_image(Q.Tissue.Image_MBR(:,:,round(Nslice/2)+1,:));

prompt = 'Pick a frame: ';
NFramePick = input(prompt);
Image_peak = Q.Tissue.Image_MBR(:,:,:,NFramePick);

for i=1:Nslice
    background_image = mean(Q.Tissue.Image_MBR(:,:,i,:),4);
    color_image = bld_flow_map(:,:,i);
    
    figure,imagesc(background_image),axis image,axis off,colormap gray
    prompt = 'Is this Normal/Apical? (0/1) ';
    flag = input(prompt);
    close
    if flag == 0
        Seg_mask{i} = Quant.segment_short_axis(color_image);
    elseif flag == 1
        Seg_mask{i} = Quant.segment_short_axis_apical(color_image);
    end
end

for i=1:Nslice
    Mask_all(:,:,i) = logical(sum(Seg_mask{i},3));
end

Q.mask.Seg_mask = Seg_mask;
Q.mask.Mask_all = Mask_all;
Q.mask.Image_show = Image_peak;

close all
clc

end


