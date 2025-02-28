function [dic_real,SI,PD_SI] = get_MISS_dic_with_phase_interleaving_2_SMS_3(flip_angle,TR,PDlines)
%--------------------------------------------------------------------------
%   [dic_real,SI,PD_SI] = get_MISS_dic_with_phase_interleaving_2_SMS_3(flip_angle,TR,PDlines)
%--------------------------------------------------------------------------
%   This function generates the dictionary used to convert signal intensity
%   time curves to [Gd] curves for blood flow quantification.
%--------------------------------------------------------------------------
%   Inputs:      
%       - flip_angle: acquisition flip angle of perfusion frames [scalar]
%       - TR: acquisition repetition time of perfusion frames [scalar]
%       - PDlines: number of proton density lines acquired [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - dic_real: dictionary for [Gd] estimation (SI/PD_SI) [nset*nsl,nT1]
%       - SI: raw dictionary signal intensity [nset*nsl,nT1]
%       - PD_SI: raw proton density dictionary signal intensity [nset*nsl,nT1]
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%           - nT1: T1 values
%--------------------------------------------------------------------------

flip_angle_PD = 2; dic_real = zeros(6,3000);

load('../mfiles/+Quant/SliceProfile/miss_excitation.mat')
slice_profile_excitation = Slice.Profile(1:1000);
slice_phase_excitation = exp(-1i.*Slice.Phase(1:1000));

load('../mfiles/+Quant/SliceProfile/miss_transition.mat')
slice_profile_transition = Slice.Profile(1:1000);

flip_angle_slice_excitation = flip_angle.*slice_profile_excitation;
flip_angle_slice_excitation_PD = flip_angle_PD.*slice_profile_excitation;

flip_angle_slice_transition = flip_angle.*slice_profile_transition;
flip_angle_slice_transition_PD = flip_angle_PD.*slice_profile_transition;

Mz = ones(4000,1); M0 = 1;

parfor Simu_T1 = 1:3000

    Slab_T1 = ones(4000,1)*Simu_T1;

    SI_one_slice = Quant.simu_steady_state_MISS_with_phase_interleaving_2_SMS_3(flip_angle_slice_excitation,flip_angle_slice_transition,slice_phase_excitation,Mz,M0,TR,Slab_T1,3000);
    SI_one_slice_PD = Quant.simu_steady_state_MISS_with_phase_interleaving_2_SMS_3(flip_angle_slice_excitation_PD,flip_angle_slice_transition_PD,slice_phase_excitation,Mz,M0,TR,Slab_T1,PDlines);

    SI(:,Simu_T1) =  SI_one_slice(:,end);
    PD_SI(:,Simu_T1) = mean(SI_one_slice_PD,2);
    
    dic_real(:,Simu_T1) = SI(:,Simu_T1) ./ PD_SI(:,Simu_T1);
end
