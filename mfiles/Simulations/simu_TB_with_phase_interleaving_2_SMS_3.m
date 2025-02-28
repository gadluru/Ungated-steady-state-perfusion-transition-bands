function [SI_one_slice,Mz] = simu_TB_with_phase_interleaving_2_SMS_3(flip_angle_slice_excitation,flip_angle_slice_transition,slice_phase_excitation,Mz,M0,TR,Slab_T1,Nalpha)
%--------------------------------------------------------------------------
%   [SI_one_slice,Mz] = simu_TB_with_phase_interleaving_2_SMS_3_with_motion(flip_angle_slice_excitation,flip_angle_slice_transition,slice_phase_excitation,Mz,M0,TR,Slab_T1,Nalpha,max_disp,heart_rate)
%--------------------------------------------------------------------------
%   This function models the Continuous Radial Interleaved
%   simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP)
%   with transition bands acquisition to generate a signal dictionary
%   simulations with no motion.
%--------------------------------------------------------------------------
%   Inputs:      
%       - flip_angle_slice_excitation: excitation slice profile with flip angle applied [nex,1]
%           - nex: length of excitation slice profile
%       - flip_angle_slice_transition: transition slice profile with flip angle applied [ntb,1]
%           - ntb: length of transition slice profile
%       - slice_phase_excitation: phase of the excitation slice profile [nex,1]
%           - nex: length of excitation slice profile
%       - Mz: magnetization history [nto,1] 
%           - nto: length of simulated acquisition
%       - M0: equilibrium magnetization [scalar]
%       - TR: acquisition repetition time of perfusion frames [scalar]
%       - Slab_T1: Slab of T1 used for dictionary estimation [nto,1]
%           - nto: length of simulated acquisition
%       - Nalpha: number of excitations [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - SI_one_slice: raw dictionary signal intensity [nset*nsl,nsa]
%           - nset: number of slice groups
%           - nsl: number of sms slices (SMS, multiband=3)
%           - nsa: number of excitations per slice
%       - Mz: magnetization history [nto,Nalpha*2+1]
%           - nto: length of simulated acquisition
%           - even indices: logitudinal magnetization after excitation
%           - odd indices: longitudinal magnetization after recovery
%--------------------------------------------------------------------------

SMS_slices = [1:1000,401:1400,801:1800] + 1000; slice_gap = 200;
transition_position = [651:1650;2351:3350];

set = repmat([1,2],[1,Nalpha/2]);
phase_excitation = vec(repmat([0,1,2],[2,Nalpha/6]))';
SI_set = zeros(2,Nalpha/2);

for i=1:Nalpha
    flip_angle_all = zeros(4000,1); slice_position = SMS_slices;

    if set(i) == 2
        slice_position = slice_position + slice_gap;
        flip_angle_all(transition_position(1,:)) = flip_angle_all(transition_position(1,:)) + flip_angle_slice_transition;
    else
        flip_angle_all(transition_position(2,:)) = flip_angle_all(transition_position(2,:)) + flip_angle_slice_transition;
    end
    
    flip_angle_all(slice_position(1:1000)) = flip_angle_all(slice_position(1:1000)) + flip_angle_slice_excitation;
    flip_angle_all(slice_position(1001:2000)) = flip_angle_all(slice_position(1001:2000)) + flip_angle_slice_excitation;
    flip_angle_all(slice_position(2001:3000)) = flip_angle_all(slice_position(2001:3000)) + flip_angle_slice_excitation;
    
    Mz(:,end+1) = Mz(:,end).*cos(flip_angle_all/180*pi);
    
    switch phase_excitation(i)
        case 0
            n_phase = [0,0,0];
        case 1
            n_phase = [0,2,1];
        case 2
            n_phase = [0,1,2];
    end
    
    N = ceil(i/2);
    for iii = [1,3,5]-1+set(i)
        n_SMS_slice = ceil(iii/2);
        
        slice_location = (1:1000) + 1000 + slice_gap*(iii-1);
        
        SI_set(set(i),N) = SI_set(set(i),N) + sum(Mz(slice_location,end-1).*sin(flip_angle_slice_excitation/180*pi).*slice_phase_excitation) .* exp(1i*n_phase(n_SMS_slice)*2*pi/3);
    end

    Mz(:,end+1) = M0 - (M0 - Mz(:,end)).*exp(-TR./Slab_T1);
    
end

SI_one_slice(1,:) = abs(SI_set(1,1:3:end) + SI_set(1,2:3:end) + SI_set(1,3:3:end));
SI_one_slice(2,:) = abs(SI_set(2,1:3:end) + SI_set(2,2:3:end) + SI_set(2,3:3:end));

SI_one_slice(3,:) = abs(SI_set(1,1:3:end) + SI_set(1,2:3:end)*exp(-1i*4*pi/3) + SI_set(1,3:3:end)*exp(-1i*2*pi/3));
SI_one_slice(4,:) = abs(SI_set(2,1:3:end) + SI_set(2,2:3:end)*exp(-1i*4*pi/3) + SI_set(2,3:3:end)*exp(-1i*2*pi/3));

SI_one_slice(5,:) = abs(SI_set(1,1:3:end) + SI_set(1,2:3:end)*exp(-1i*2*pi/3) + SI_set(1,3:3:end)*exp(-1i*4*pi/3));
SI_one_slice(6,:) = abs(SI_set(2,1:3:end) + SI_set(2,2:3:end)*exp(-1i*2*pi/3) + SI_set(2,3:3:end)*exp(-1i*4*pi/3));

end