function step_size = line_search_ADMM(old,update,Data,para)
%--------------------------------------------------------------------------
%   step_size = line_search_ADMM(old, update, Data, para)
%--------------------------------------------------------------------------
%   This function trys to find a suitable step size to perform a CG update 
%   for ADMM reconstruction. The function starts with a step size adopted 
%   from last iteration, and multiply it by 1.3 (magic number). If the step
%   size yeilds a cost that is larger than the previous cost, it shrinks 
%   the step size by 0.8 (magic number again). If it yeilds a cost that is 
%   smaller than the previous cost, it will increase the step size by 1.3 
%   until it no longer yeild a smaller cost. The maximum number of trys is 15.
%--------------------------------------------------------------------------
%   Inputs:
%       - old: image series from previous iteration [nx,ny,nt,nset*nsl]
%       - update: conjugate gradient update [nx,ny,nt,nset*nsl]
%       - Data: Reconstruction variables after cardiac phase fixing [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,nset,nsl]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt,1,1,nset,nsl]
%           - SMS: CAIPI phase modulation matrix [1,1,1,1,nsl,1,nsl]
%           - first_guess: preliminary STCR [nx,ny,nt,nset*nsl]
%           - Y: ADMM Lagrangian multiplier (following first iteration of ADMM Step 1) [nx,ny,nt,nset*nsl]
%           - sens_map: coil sensitivity maps [nx,ny,1,nc,nset,nsl]
%           - ramp_filter: radial sampling filter [nx,ny,nt,1,1,nset]
%           - sms_filter: SMS acquisition filter [nx,ny,nt,1,1,nset]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames after cardiac phase fixing
%               - nc: number of PCA coils
%               - nset: number of slice groups
%               - nsl: number of sms slices (SMS, multiband=3)
%           - llr: Motion tracked patches
%               - Npatch: number of motion tracked patches
%               - idx: linear indices of motion tracked patches [prod(patch_size),Ncycle,Npatch]
%               - mask_intensity: averaging mask of motion tracked patches [nx,ny,nset*nsl,nt]
%               - mask: location mask of motion tracked patches [nx,ny,nset*nsl,nt]
%                   - nx: number of measurements in spatial x-dimension
%                   - ny: number of measurements in spatial y-dimension
%                   - nt: number of time frames after cardiac phase fixing
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%                   - Ncycle: number of cardiac cycles in the image series
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - step_size: conjugate gradient step size to compute updated image
%                    series for ADMM reconstruction [scalar]
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

step_start = para.Recon.step_size(end)*1.3;% magic number
tau = 0.8; max_try = 15; step_size = step_start; flag = 0;

cost_old = para.Cost.totalCost(end);

for i=1:max_try

    new = old + step_size*update;
    fidelity_norm = compute_fidelity_for_line_search_yt(sr_forward(new,Data.size),Data,para);
    
    if isfield(Data, 'Y')
        [~,cost_new] = Cost_STCR_step_3(new, Data.first_guess + Data.Y, fidelity_norm, Data.llr.MT, para);
    else
        [~,cost_new] = Cost_STCR_step_2(new, fidelity_norm, Data.llr.MT, para);
    end
    
    if cost_new > cost_old && flag == 0
        step_size = step_size*tau;
    elseif cost_new < cost_old 
        step_size = step_size*1.3;
        cost_old = cost_new;
        flag = 1;
    elseif cost_new > cost_old && flag == 1
        step_size = step_size/1.3;
        return
    end
    
end
end