function Q = get_MISS_tissue_2_SMS_3(Q)
%--------------------------------------------------------------------------
%   Q = get_MISS_tissue_2_SMS_3(Q)
%--------------------------------------------------------------------------
%   Quantification for ungated, free-breathing, cardiac phase-resolved
%   myocardial perfusion MRI using Continuous Radial Interleaved
%   simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP). 
%   This quantification code builds upon the CRIMP framework
%   to perform quantification for the same acquisition with the addition of
%   transition bands to improve cardiac perfusion. This code uses the pTV
%   registration toolbox for motion compensation and Bullseye Plot code
%   from the Matlab File Exchange.
%--------------------------------------------------------------------------
%   Inputs:      
%       - Q: Quantification parameters [struct]
%           -dir: file location of reconstruction for quantification
%--------------------------------------------------------------------------
%   Outputs:
%       - Q: Quantification parameters [struct]
%           -dir: file location of reconstruction for quantification
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
%           -AIF
%               -mask: binary mask in LV blood pool [nx,ny,nset*nsl]
%               -SI_curve_second_all: raw signal AIF curve [nset*nsl,nt]
%               -T1_curve_second_all: T1 AIF curve [nset*nsl,nt]
%               -Gd_curve_second_all: [Gd] AIF curve [nset*nsl,nt]
%               -slice_pick: slices used for AIF (averaged if more than one slice) [scalar/array]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%           -param
%               -TD: time spacing between perfusion frames after interpolation [scalar]
%               -TimeStamps: time duration of perfusion after interpolation [1,nt]
%               -PDlines: number of proton density lines acquired [scalar]
%               -FlipAngle: acquisition flip angle of perfusion frames
%               -TR: acquisition repetition time of perfusion frames
%               -Nprecontrast: number of precontrast frames after interpolation
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
%               -bld_flow: Mean/SD AHA regional blood flow values [6,6,2]
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
%   Reference:
%       [1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G.
%           Improving ungated steady-state cardiac perfusion using transition
%           bands. Magn Reson Med 2025.
%       [2] Tian Y, Mendes J, Wilson B, Ross A, Ranjan R, DiBella E, Adluru G.
%           Whole-heart, ungated, free-breathing, cardiac-phase-resolved 
%           myocardial perfusion MRI by using Continuous Radial Interleaved
%           simultaneous Multi-slice acquisitions at sPoiled steady-state 
%           (CRIMP). Magn Reson Med 2020;84(6):3071-3087
%       [3] Vishnevskiy V, Gass T, Szekely G, Tanner C, Goksel O. Isotropic
%           total variation regularization of displacements in parametric image
%           registration. IEEE Trans Med Imaging 2017;36(2):385-395
%       [4] Adrian (2025). Bullseye Plot.zip (https://www.mathworks.com/
%           matlabcentral/fileexchange/47454-bullseye-plot-zip),
%           MATLAB Central File Exchange. Retrieved February 26, 2025. 
%--------------------------------------------------------------------------
%   Author:
%       Johnathan Le
%       E-mail: le.johnv@outlook.com
%--------------------------------------------------------------------------

load(fullfile(Q.dir.file.folder,Q.dir.file.name),'Image_sys','Image_PD','para');

Image = abs(Image_sys);
Image_PD = abs(Image_PD);

[nx,ny,nt,nsl] = size(Image);

%% Parameter Initialization
siz = size(Image_PD);
if length(siz)>3
    if siz(4)<siz(3)
        Image_PD = permute(Image_PD,[1,2,4,3]);
    end
end

siz = size(Image);
if length(siz)>3
    if siz(4)<siz(3)
        Image = permute(Image,[1,2,4,3]);
    end
end

Image = zero_pad_image(Image);
Image_PD = zero_pad_image(Image_PD);
motion_mask = zero_pad_image(para.Recon.self_gating.motion_mask);
cardiac_mask = zero_pad_image(para.Recon.self_gating.cardiac_mask);

for i=1:nsl
    motion_mask(:,:,i) = imfill(motion_mask(:,:,i) > 0.9,'holes');
    cardiac_mask(:,:,i) = imfill(cardiac_mask(:,:,i) > 0.9,'holes');
end
motion_mask = logical(squeeze(motion_mask)); cardiac_mask = logical(squeeze(cardiac_mask));

center = floor(regionprops(bwconncomp(any(cardiac_mask,3)),'Centroid').Centroid);

%% Perfusion Registration

fprintf('Perfusion Registration...\n')

[nx,ny,nsl,nt] = size(Image);

tic
Image_MBR = zeros(size(Image));
for i=1:nsl
    Image_MBR(:,:,i,:) = Quant.pTV_perf(squeeze(Image(:,:,i,:)));
end

for i=1:nsl
    Image_MBR(:,:,i,:) = Quant.mbr(squeeze(Image_MBR(:,:,i,:)),0);
end
toc

%% Proton Density Registration

fprintf('Proton Density Registration...\n')

tic
Image_deform_MC_PD = zeros(nx,ny,nsl);
for i=1:nsl
    Image_deform_MC_PD(:,:,i) = Quant.pTV_PD_reg(squeeze(Image_PD(:,:,i,:)));
end
toc

%% check registration
line_profile_x = squeeze(cat(4,Image_deform_MC_PD(center(2),:,:,:),Image_MBR(center(2),:,:,:)));
line_profile_x = permute(line_profile_x,[1,3,2]);

line_profile_y = squeeze(cat(4,Image_deform_MC_PD(:,center(1),:,:),Image_MBR(:,center(1),:,:)));
line_profile_y = permute(line_profile_y,[1,3,2]);

f1 = figure; figure(f1); movegui('center')
for i=2:nsl-1
    subplot(2,nsl-2,i-1)
    imagesc(line_profile_x(:,:,i))
    subplot(2,nsl-2,i+nsl-3)
    imagesc(line_profile_y(:,:,i))
end

Q.Tissue.Image_PD = mean(Image_deform_MC_PD,4);
Q.Tissue.Image_PD(Q.Tissue.Image_PD<mean(Q.Tissue.Image_PD(:))/10) = mean(Q.Tissue.Image_PD(:))/10;

Q.Tissue.Image_MBR = Image_MBR(:,:,:,3:end-2);

[nx,ny,nsl,nt] = size(Q.Tissue.Image_MBR);

%% get AIF mask

if ismember('AIFmask',who('-file','mask.mat'))
    Q.AIF.mask = load('mask.mat','AIFmask').AIFmask;
else
    fprintf('Draw ROI of the LV to calculate AIF...\n');

    [~,frame] = max(squeeze(sum(sum(Q.Tissue.Image_MBR .* cardiac_mask))),[],2);

    f2 = figure; figure(f2); movegui('center')
    for i=1:nsl
        imagesc(mean(Q.Tissue.Image_MBR(:,:,i,frame(i)-1:frame(i)+1),4)),axis off,axis equal,colormap gray
        mask_all(:,:,i,1) = roipoly;
    end

    for i=1:nsl
        idx = Q.Tissue.Image_MBR(:,:,i,frame(i)) .* mask_all(:,:,i,1);
        Q.AIF.mask(:,:,i,1) = idx > max(vec(idx))*0.8;
    end
    close(f2); clear f2
end

%% Time Stamp

Q.param.TD = 1;

TimeStamps = para.kSpace_info.TimeStamp(3:end-2);
TimeStamps = TimeStamps - TimeStamps(1);

Q.param.TimeStamps = 0:Q.param.TD:floor(max(TimeStamps));

[nx,ny,nsl,nt] = size(Q.Tissue.Image_MBR);

Q.Tissue.Image_MBR = interp1(TimeStamps,reshape(Q.Tissue.Image_MBR,[nx*ny*nsl,nt])',Q.param.TimeStamps);
Q.Tissue.Image_MBR = reshape(Q.Tissue.Image_MBR',[nx,ny,nsl,size(Q.Tissue.Image_MBR,1)]);

%% load parameter

fprintf('Estimating Dictionary...\n')

Q.param.PDlines = para.kSpace_info.NumberOfPDlines;
Q.param.FlipAngle = para.kSpace_info.Protocol.adFlipAngleDegree;
Q.param.TR = para.kSpace_info.Protocol.alTR(2)/1000;

tic
[dic,SI,PD_SI] = Quant.get_MISS_dic_with_phase_interleaving_2_SMS_3(Q.param.FlipAngle,Q.param.TR,Q.param.PDlines);
toc

%% AIF GD curve

fprintf('AIF Signal Intensity to Gadolinium...\n')

SI_curve = squeeze(sum(sum(Q.Tissue.Image_MBR.*Q.AIF.mask)));
SI_curve = SI_curve./squeeze(sum(sum(Q.AIF.mask)));

PD_SI_curve = squeeze(sum(sum(Q.Tissue.Image_PD.*Q.AIF.mask)));
PD_SI_curve = PD_SI_curve./squeeze(sum(sum(Q.AIF.mask)));

Q.param.Nprecontrast = Quant.auto_firstpass(mean(SI_curve,1)); Q.param.Nprecontrast = 1:Q.param.Nprecontrast;

SI_ratio = SI_curve ./ PD_SI_curve; SI_ratio_pre = mean(SI_ratio(:,Q.param.Nprecontrast),2);

tic
for i=1:nsl
    [Gd_curve(i,:),T1_curve(i,:)] = Quant.SI2Gd_dic_real(SI_ratio(i,:),SI_ratio_pre(i,:),dic(i,:),SI(i,:),PD_SI(i,:),Q.param.Nprecontrast);
end
toc

Q.AIF.SI_curve_second_all = SI_curve;
Q.AIF.Gd_curve_second_all = Gd_curve;
Q.AIF.T1_curve_second_all = T1_curve;

Quant.plotAIF(Q.AIF.Gd_curve_second_all',Q.param.TimeStamps)

Q.AIF.slice_pick = [];
while ~isnumeric(Q.AIF.slice_pick) || isempty(Q.AIF.slice_pick)
    Q.AIF.slice_pick = 5;%input('Select AIF slices: ');
end

clear SI_curve PD_SI_curve SI_ratio SI_ratio_pre dSI_curve

%% Tissue GD curve

fprintf('Tissue Signal Intensity to Gadolinium...\n')

SI_ratio = Q.Tissue.Image_MBR ./ Q.Tissue.Image_PD; SI_ratio_pre = mean(SI_ratio(:,:,:,Q.param.Nprecontrast),4);

tic
for i=1:nsl
    for j=1:nx
        [Gd_all(j,:,i,:),T1_all(j,:,i,:)] = Quant.SI2Gd_dic_tissue_real(squeeze(SI_ratio(j,:,i,:)),squeeze(SI_ratio_pre(j,:,i,:)),dic(i,:),SI(i,:),PD_SI(i,:),Q.param.Nprecontrast);
    end
end
toc

Q.Tissue.SI_curve_second_all = Q.Tissue.Image_MBR;
Q.Tissue.T1_curve_second_all = T1_all;
Q.Tissue.Gd_curve_second_all = Gd_all;

%% model fit

fprintf('Blood Flow Quantification...\n')

nsecond = size(Q.Tissue.Gd_curve_second_all,4);

Q.Tissue.Mask_heart = motion_mask;

Tissue_Gd_curve_second_all = Q.Tissue.Gd_curve_second_all(repmat(Q.Tissue.Mask_heart,[1,1,1,nsecond]));
Tissue_Gd_curve_second_all = reshape(Tissue_Gd_curve_second_all,[sum(Q.Tissue.Mask_heart(:)),nsecond])';

AIF_Gd_curve_second_all = mean(Q.AIF.Gd_curve_second_all(Q.AIF.slice_pick,:),1);

idx_throw = Tissue_Gd_curve_second_all > 2 | Tissue_Gd_curve_second_all < -0.2;
idx_throw = logical(sum(idx_throw,1));

idx = find(Q.Tissue.Mask_heart);
Q.Tissue.Mask_heart(idx(idx_throw)) = false;
Tissue_Gd_curve_second_all(:,idx_throw) = [];

tic
Q = Quant.get_2CM(AIF_Gd_curve_second_all,Tissue_Gd_curve_second_all,Q);
toc

%%

if ismember('mask',who('-file','mask.mat'))
    Q.mask = load('mask.mat','mask').mask;
else
    Q = Quant.get_mask_AHA(Q);
end

Q = Bullseye(Q);

save('Quant.mat','Q','-v7.3')

return

