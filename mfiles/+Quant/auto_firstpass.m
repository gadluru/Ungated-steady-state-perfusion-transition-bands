function [precontrast,firstpass_cutoff]=auto_firstpass(curve)
%--------------------------------------------------------------------------
%   [precontrast,firstpass_cutoff]=auto_firstpass(curve)
%--------------------------------------------------------------------------
%   This function returns the index of the last precontrast time frame and
%   the index of last time frame of the first pass of contrast.
%--------------------------------------------------------------------------
%   Inputs:      
%       - curve: arterial input function [1,nt]
%           -nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - precontrast: index of last precontrast time frame [scalar]
%       - firstpass_cutoff: index of last time frame of the first pass of
%       contrast
%--------------------------------------------------------------------------

h = fspecial('gaussian', [1 size(curve,2)], 4);
filt_aif = imfilter(curve,h,'replicate');

diff_aif=diff(filt_aif);
step_var=diff_aif;
step_var(find(step_var<0))=-1;
step_var(find(step_var>=0))=1;
precontrast = find(diff_aif>0.05);
precontrast = precontrast(1);

min_diff=find(diff_aif==min(diff_aif));
tmp1=find(step_var>0);

try
tmp2=tmp1(find(tmp1>=min_diff));
firstpass_cutoff=tmp2(1)-1;
catch
firstpass_cutoff=length(curve);    
end
