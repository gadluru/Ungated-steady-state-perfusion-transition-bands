
%% Improving Ungated Steady-State Perfusion using Transition Bands
% This code performs the reconstruction from the MRM publication,
% Mendes JK, Le JV, Arai AE, et al. Magn Reson Med 2025. This code
% was developed and tested on a RockLinux 8.6 operating system, with an AMD
% EPYC Milan 7543 32 2.8GHz 256 MB cache, 512 GB RAM and Nvidia A100 gpus.
% Systems with ~64 GB RAM may encounter issues with insufficient memory 
% that require code changes to reduce memory requirements during certain 
% function calls.
%--------------------------------------------------------------------------
%   Reference:
%       [1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G.
%           Improving ungated steady-state cardiac perfusion using transition
%           bands. Magn Reson Med 2025.
%--------------------------------------------------------------------------
%   Author:
%       Johnathan Le
%       E-mail: le.johnv@outlook.com
%--------------------------------------------------------------------------

addpath(genpath('mfiles/'))

if ~exist('RawData','dir')
    mkdir('RawData')
end

% example human dataset download
if ~isfile('RawData/example_human_dataset.mat')
    url = 'https://dataverse.harvard.edu/api/access/datafile/10907539';
    websave('RawData/example_human_dataset.mat',url);
end

all_mat = dir('RawData/example_human_dataset.mat');

for i=1
    %% set reconstruction parameters

    para.dir.load_kSpace_name = all_mat(i).name;
    para.dir.load_kSpace_dir = [all_mat(i).folder,'/'];
    
    load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name])
    
    [nx,nor_all,nc] = size(kSpace);

    para.setting.ifGPU = 1; % This code requires ~50 Gb of GPU memory, set flag to 0 to run on CPU

    para.Recon.noi = 50; % number of iterations for reconstruction. This is changed several times during the reconstruction.
    para.Recon.kSpace_center = floor(nx/2) + 1; % kspace center of readout

    para.Recon.nor = 12; % Number of rays per frame, temporal footprint, should be multiple of 6 to account for SMS/Transition Band phase modulation.
    para.Recon.nor_sl = 12; % Number of rays between ecah frame, temporal resolution, should be multiple of 6 to account for SMS/Transition Band phase modulation.
    
    para.weight_sTV = 0.00; % spatial TV constraint regularization weight
    para.weight_slTV = 0.00; % slice TV constraint regularization weight
    para.weight_tTV = 0.025; % temporal TV constraint regularization weight   

    para.Recon.ramp_filter = 1; % radial sampling filter
    para.Recon.sms_filter = 1; % SMS acquisition filter

    para.Recon.weight_l2 = 0.01; % LLR constraint regularization weight for ADMM Perfusion reconstruction
    para.Recon.weight_lr = 0.01; % LLR constraint regularizaiton weight for Proton Density reconstruction
    
    para.Recon.bloc_x = 8; % X-direction patch length for LLR constraint
    para.Recon.bloc_y = 8; % Y-direction patch length for LLR constraint
    para.Recon.tau = 20; % LLR parameter for soft thresholding
    
    para = prepare_para(para);
    
    %% run the main reconstruction function
    recon_2D_ungated_SPGR_nSMS_ADMM(kSpace,para);
end
