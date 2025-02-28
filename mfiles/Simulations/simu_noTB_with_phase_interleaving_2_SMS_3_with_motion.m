function SI_one_slice = simu_noTB_with_phase_interleaving_2_SMS_3_with_motion(flip_angle_slice_excitation,slice_phase_excitation,Mz,M0,TR,Slab_T1,Nalpha,max_disp,heart_rate)
%--------------------------------------------------------------------------
%   SI_one_slice = simu_noTB_with_phase_interleaving_2_SMS_3_with_motion(flip_angle_slice_excitation,slice_phase_excitation,Mz,M0,TR,Slab_T1,Nalpha,max_disp,heart_rate)
%--------------------------------------------------------------------------
%   This function models the Continuous Radial Interleaved
%   simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP)
%   with no transition bands acquisition to generate a signal dictionary for
%   motion simulations
%--------------------------------------------------------------------------
%   Inputs:      
%       - flip_angle_slice_excitation: excitation slice profile with flip angle applied [nex,1]
%           - nex: length of excitation slice profile
%       - slice_phase_excitation: phase of the excitation slice profile [nex,1]
%           - nex: length of excitation slice profile
%       - Mz: magnetization history [nto,1] 
%           - nto: length of simulated acquisition
%       - M0: equilibrium magnetization [scalar]
%       - TR: acquisition repetition time of perfusion frames [scalar]
%       - Slab_T1: Slab of T1 used for dictionary estimation [nto,1]
%           - nto: length of simulated acquisition
%       - Nalpha: number of excitations [scalar]
%       - max_disp: maximum slice displacement from motion simulation in units of slice thickness [scalar]
%       - heart_rate: simulated heart rate for sinusoidal motion model [scalar]
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

SMS_slices = [1:1000,401:1400,801:1800] + 400; slice_gap = 200;

T = 1:Nalpha;
total_time = T(end)*TR/1000;
Ncycle = total_time/60*heart_rate;
max_disp = 200*max_disp;
Motion = round(sin(T/T(end)*Ncycle*2*pi)*max_disp);

set = repmat([1,2],[1,Nalpha/2]);
phase = vec(repmat([0,1,2],[2,Nalpha/6]))';
SI_set = zeros(2,Nalpha/2);

for i=1:Nalpha
    flip_angle_all = zeros(3000,1); slice_position = SMS_slices;

    if set(i) == 2
        slice_position = slice_position + slice_gap;
    end
    slice_position = slice_position + Motion(i);
    
    flip_angle_all(slice_position(1:1000)) = flip_angle_all(slice_position(1:1000)) + flip_angle_slice_excitation;
    flip_angle_all(slice_position(1001:2000)) = flip_angle_all(slice_position(1001:2000)) + flip_angle_slice_excitation;
    flip_angle_all(slice_position(2001:3000)) = flip_angle_all(slice_position(2001:3000)) + flip_angle_slice_excitation;

    Mz(:,end+1) = Mz(:,end).*cos(flip_angle_all/180*pi);
    
    switch phase(i)
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

        slice_location = (1:1000) + slice_gap*(iii-1) + Motion(i) + 400;

        SI_set(set(i),N) = SI_set(set(i),N) + sum(Mz(slice_location,end-1).*sin(flip_angle_slice_excitation/180*pi).*slice_phase_excitation) .*exp(1i*n_phase(n_SMS_slice)*2*pi/3);
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