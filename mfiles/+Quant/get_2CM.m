
function Q = get_2CM(AIF,Tissue,Q)
%--------------------------------------------------------------------------
%   Q = get_2CM(AIF,Tissue,Q)
%--------------------------------------------------------------------------
%   This function estimates the model parameters from the two-compartment
%   model 
%--------------------------------------------------------------------------
%   Inputs:      
%       - AIF: arterial input function for blood flow quantification [1,nt]
%           - nt: number of perfusion time frames
%       - Tissue: tissue curves for blood flow quantification [nt,nb]
%           - nt: number of perfusion time frames
%           - nb: number of tissue curves
%       - Q
%           - Tissue
%               -Mask_heart: Mask for blood flow quantification [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%           -param
%               -TD: time spacing between perfusion frames after interpolation [scalar]
%               -TimeStamps: time duration of perfusion after interpolation [1,nt]
%                   - nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - Q
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
%--------------------------------------------------------------------------
%   Reference:
%       [1] Tofts PS, Brix G, Buckley DL, Evelhoch JL, Henderson E, Knopp
%       MV, Larsson HB, Lee TY, Mayr NA, Parker GJ, Port RE, Taylor J,
%       Weisskoff RM. Estimating kinetic parameters from dynamic
%       contrast-enhanced T1-weighted MRI of a diffusable tracer:
%       Standardized quantities and symbols. J Magn Reson Imaging
%       1999;10(3):223-232.
%--------------------------------------------------------------------------

    Q.TCM.para_low = [0.1,  0.1,  0.01,  0.01];
    Q.TCM.para_up =  [10,   30,   5,     0.5];
    Q.TCM.para_init =[0.8,  4,    0.05,  0.05];

    param_all = Quant.Fit_2CM(AIF,Tissue,Q);

    Q.TCM.bld_flow_map = single(zeros(size(Q.Tissue.Mask_heart)));
    Q.TCM.bld_flow_map(Q.Tissue.Mask_heart) = param_all(1,:);

    Q.TCM.kep_map = single(zeros(size(Q.Tissue.Mask_heart)));
    Q.TCM.kep_map(Q.Tissue.Mask_heart) = param_all(2,:);

    Q.TCM.vb_map = single(zeros(size(Q.Tissue.Mask_heart)));
    Q.TCM.vb_map(Q.Tissue.Mask_heart) = param_all(3,:);

    Q.TCM.dt_map = single(zeros(size(Q.Tissue.Mask_heart)));
    Q.TCM.dt_map(Q.Tissue.Mask_heart) = param_all(4,:);

end