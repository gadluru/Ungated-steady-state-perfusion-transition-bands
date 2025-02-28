function SI = get_TB_dic_with_phase_interleaving_2_SMS_3_with_motion(T1,flip_angle,TR,Nalpha,max_disp,heart_rate)
%--------------------------------------------------------------------------
%   SI = get_TB_dic_with_phase_interleaving_2_SMS_3_with_motion(T1,flip_angle,TR,max_disp,heart_rate)
%--------------------------------------------------------------------------
%   This function generates the signal dictionary used to simulate the Continuous 
%   Radial Interleaved simultaneous Multi-slice acquisition at sPoiled
%   steady-state (CRIMP) with transition bands using a sinusoidal motion
%   model.
%--------------------------------------------------------------------------
%   Inputs:      
%       - flip_angle: acquisition flip angle of perfusion frames [scalar]
%       - TR: acquisition repetition time of perfusion frames [scalar]
%       - Nalpha: number of excitations [scalar]
%       - max_disp: maximum slice displacement from motion simulation in units of slice thickness [scalar]
%       - heart_rate: simulated heart rate for sinusoidal motion model [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - SI: raw dictionary signal intensity [nset*nsl,nT1]
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%           - nT1: T1 values
%--------------------------------------------------------------------------

load('../../mfiles/+Quant/SliceProfile/miss_excitation.mat')
slice_profile_excitation = Slice.Profile(1:1000);
slice_phase_excitation = exp(-1i.*Slice.Phase(1:1000));

load('../../mfiles/+Quant/SliceProfile/miss_transition.mat')
slice_profile_transition = Slice.Profile(1:1000);

flip_angle_slice_excitation = flip_angle.*slice_profile_excitation;
flip_angle_slice_transition = flip_angle.*slice_profile_transition;

Mz = ones(4000,1); M0 = 1; SI = zeros(6,Nalpha/6);
parfor i=1:length(T1)

    Slab_T1 = ones(4000,1)*T1(i);

    SI(:,:,i) = simu_TB_with_phase_interleaving_2_SMS_3_with_motion(flip_angle_slice_excitation,flip_angle_slice_transition,slice_phase_excitation,Mz,M0,TR,Slab_T1,Nalpha,max_disp,heart_rate);

end
