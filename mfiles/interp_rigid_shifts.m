function shifts_interp = interp_rigid_shifts(shifts,respiratory_signal,reference_points)
%--------------------------------------------------------------------------
%   shifts_interp = interp_rigid_shifts(shifts,respiratory_signal,reference_points)
%--------------------------------------------------------------------------
%   Function to linearly fit the diastolic rigid translations to the respiration
%   signal to obtain interpolated translations to other cardiac phases.
%--------------------------------------------------------------------------
%   Inputs:      
%       - shifts: calculated x-, y-, and z- rigid shifts [3,noi,nt]
%           - noi: number of iterations
%           - nt: number of diastolic time frames
%       - respiratory_signal: extracted respiration signal [nt,1]
%           - nt: number of time frames
%       - reference_points: location of diastolic time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - shifts_interp: calculated x-,y-, and z- rigid shifts [3,nt]
%           - nt: total number of time frames
%--------------------------------------------------------------------------

shifts_interp = zeros(size(shifts,1),length(respiratory_signal));
shift_ref = respiratory_signal(reference_points);

shifts = squeeze(sum(shifts,2));
for i=1:size(shifts,1)
    f = fit(double(shift_ref), double(shifts(i,:)'), 'poly1' );
    shifts_interp(i,:) = f(respiratory_signal);
end