function count_mask_voxels()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Instructions: 
% 
% 1. Run function by placing it in your MATLAB environment path and calling
% 'count_mask_voxels'
% 
% 2. When prompted, select the directory containing the root path to your
% input files (note that input subpaths likely need to change 
% to reflect your own folder structures): 
% 
% wholehead_file: NIFTI of your original non-skull-stripped volume
% brain_mask_file: NIFTI of brain mask (note that code currently accounts for an up-down flip
% in this relative to the wholehead file) 
% ventricle_mask_file: NIFTI of brain ventricle mask to be subtracted for final
% brain volume calculation 
% 
% 3. Function will visualize your skull-stripped brain with overlaid ventricle mask for sanity check. 
% Output file 'xxx_MaskVoxelCounts.csv' will give slicewise voxel counts for 
% skull-stripped brain and ventricle masks that can be used in final brain volume
% calculation 
%
% Author: Kelley Swanberg (Lunds universitet, 2024) 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % count_mask_voxels 
    root_folder = uigetdir(); 
    
    % Define input file locations 
    wholehead_file = sprintf('%s//T2_TurboRARE_highres.nii', root_folder); 
    brain_mask_file =  sprintf('%s//BrainMask//BrainMask.nii.gz', root_folder); 
    ventricle_mask_file = sprintf('%s//BrainMask//T2_TurboRARE_highres_VentricleMask_v2.nii', root_folder);  
    
    % Load input files for whole head, brain mask, and ventricle mask 
    wholehead = niftiread(wholehead_file); 
    brain_mask_flipped = niftiread(brain_mask_file); 
    brain_mask = int16(imbinarize(flipud(brain_mask_flipped))); % Need to flip relative to original head and convert to uint16
    ventricle_mask = niftiread(ventricle_mask_file); 
    brain = brain_mask .* wholehead; 
    vol_info = niftiinfo(wholehead_file); 
    
    % Output CSV for slicewise voxel counts
    combined_table = []; 
    var_names = {'Slice_Number', 'Brain_Voxel_Num', 'Ventricle_Voxel_Num'};
    
    numslices = min(size(wholehead)); 
    
    output_image_directory = 'Mask_Images'; 
    output_image_path = sprintf('%s//%s', root_folder, output_image_directory); 
    mkdir(output_image_path); 
    
    coeff = 1.5 * max(max(max(wholehead))); 
    
    for ii=1:numslices
        slice_number = ii; 
        brain_voxels = sum(sum(brain_mask(:,:,ii))); 
        ventricle_voxels = sum(sum(ventricle_mask(:,:,ii))); 
    
        combined_table(ii, 1) = slice_number;
        combined_table(ii, 2) = brain_voxels;
        combined_table(ii, 3) = ventricle_voxels;
    
        image_to_plot = wholehead(:, :, ii) + coeff * int16(edge(brain_mask(:, :, ii))) + coeff * int16(edge(ventricle_mask(:, :, ii))); 
        colormap(gray); 
        imagesc(rot90(image_to_plot)); 
    
        png_name = sprintf('%s//%s//Slice_%d.png', root_folder, output_image_directory, ii); 
        saveas(gcf, png_name, 'png'); 
    end 
    
    % Add row names 
    combined_table = array2table(combined_table); 
    combined_table.Properties.VariableNames = var_names; 
    
    % Write final table to combined CSV
    [case_unclean, ~] = fileparts(root_folder ); 
    [~, case_clean] = fileparts(case_unclean); 
    csv_filename = sprintf('%s_MaskVoxelCounts.csv', case_clean); 
    writetable(combined_table, csv_filename, 'WriteRowNames', true); 
    
    output_fullpath = sprintf('%s//%s_MaskedBrain.nii', root_folder, case_clean);
    niftiwrite(int16(brain),output_fullpath,vol_info); 

end
