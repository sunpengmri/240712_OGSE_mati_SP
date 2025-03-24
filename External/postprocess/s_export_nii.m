function s_export_nii(ref_nii,map,out_file)

% s_export_nii Summary of this function goes here
%   input: 
%        ref_nii: path and name of nii reference
%        map: the quantative map
%        out_file: path and name of map.nii 
    [size_x,size_y]=size(map);
    nii_info=niftiinfo(ref_nii); 
    nii_info.Datatype='double';
    nii_info.AdditiveOffset = 0 ; 
    nii_info.raw.scl_slope = 1;
    nii_info.raw.scl_inter = 0;
    map_nii=flipud(rot90(map));
%     if size(map_nii,3)>1
   % map_nii=fliplr(map_nii);
    niftiwrite(map_nii,out_file,nii_info);
%     else
%         niftiwrite(reshape(map_nii,[size_x size_y 1]),out_file,nii_info);
%     end
end