function [x_window_begin,x_window_end,y_window_begin,y_window_end,z_window_begin,z_window_end] = get_search_window(x_begin,y_begin,z_begin,search_size,patch_size,offset_cycle,nx,ny,nsl)
%--------------------------------------------------------------------------
%   [x_window_begin,x_window_end,y_window_begin,y_window_end,z_window_begin,z_window_end] = get_search_window(x_begin,y_begin,z_begin,search_size,patch_size,offset_cycle,nx,ny,nsl)
%--------------------------------------------------------------------------
%   This function returns the pixel locations of the search window for
%   possible locations a patch can move
%--------------------------------------------------------------------------
%   Inputs:      
%       - x_begin: X-direction starting location of patch [scalar]
%       - y_begin: Y-direction starting location of patch [scalar]
%       - z_begin: Z-direction starting location of patch [scalar]
%       - search_size: possible search size for patch tracking [9,9,4]
%       - patch_size: 3D patch size [5,5,2]
%       - offset_cycle: interpolated diastolic rigid translations for guided patch tracking [3,nt]
%       - nx: number of measurements in spatial x-dimension
%       - ny: number of measurements in spatial y-dimension
%       - nt: number of time frames after cardiac phase fixing
%       - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------
%   Outputs:
%       - x_window_begin: X-direction starting location of search window [scalar]
%       - x_window_end: X-direction ending location of search window [scalar]
%       - y_window_begin: Y-direction starting location of search window [scalar]
%       - y_window_end: Y-direction ending location of search window [scalar]
%       - z_window_begin: Z-direction starting location of search window [scalar]
%       - z_window_end: Z-direction ending location of search window [scalar]
%--------------------------------------------------------------------------

x_window_begin = max(x_begin-round((search_size(1)-patch_size(1))/2)+offset_cycle(1),1); x_window_end = min(x_window_begin-1+search_size(1),nx);
if x_window_begin==1
    x_window_end = max(x_window_end,x_window_begin+patch_size(1)-1);
end
if x_window_end==nx
    x_window_begin = min(x_window_begin,x_window_end-patch_size(1)+1);
end

y_window_begin = max(y_begin-round((search_size(2)-patch_size(2))/2)+offset_cycle(2),1); y_window_end = min(y_window_begin-1+search_size(2),ny);
if y_window_begin==1
    y_window_end = max(y_window_end,y_window_begin+patch_size(2)-1);
end
if y_window_end==ny
    y_window_begin = min(y_window_begin,y_window_end-patch_size(2)+1);
end

z_window_begin = max(z_begin-round((search_size(3)-patch_size(3))/2)+offset_cycle(3),1); z_window_end = min(z_window_begin-1+search_size(3),nsl);
if z_window_begin==1
    z_window_end = max(z_window_end,z_window_begin+patch_size(3)-1);
end
if z_window_end==nsl
    z_window_begin = min(z_window_begin,z_window_end-patch_size(3)+1);
end

end