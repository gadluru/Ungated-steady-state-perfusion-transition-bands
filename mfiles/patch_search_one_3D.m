function offset = patch_search_one_3D(patch,search_window)
%--------------------------------------------------------------------------
%   offset = patch_search_one_3D(patch,search_window)
%--------------------------------------------------------------------------
%   This function finds the location of the patch in the next cardiac cycle
%   based on the minimum difference of the patch in the search window
%--------------------------------------------------------------------------
%   Inputs:      
%       - patch: 3D patch to be tracked [5,5,2]
%       - search_window: 3D search window of possible tracking locations [9,9,4]
%--------------------------------------------------------------------------
%   Outputs:
%       - offset: 3D offset between the patch and its position in the next cardiac cycle [x,y,z]
%--------------------------------------------------------------------------

patch_size = size(patch);
search_size = size(search_window);
Npatches = search_size - patch_size + 1;

patches_all = zeros([patch_size,Npatches]);

for i=1:Npatches(1)
    for j=1:Npatches(2)
        for k=1:Npatches(3)
            patches_all(:,:,:,i,j,k) = search_window(i:i+patch_size(1)-1,j:j+patch_size(2)-1,k:k+patch_size(3)-1);
        end
    end
end

d = patch - patches_all;
d = sum(sum(sum(abs(d))));
[~,idx] = min(d(:));
[x,y,z] = ind2sub(Npatches,idx);
offset = [x,y,z];

end