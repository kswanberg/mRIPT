%function Manually_Mask_ROI()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Instructions: 
% 
% 1. Run function by placing it in your MATLAB environment path and calling
% 'Manually_Mask_ROI'
% 
% 2. When prompted, select the Nifti file on which you would like to draw your
% slicewise ROIs
% 
% 3. Draw up to 'roi_num_per_slice' ROI shapes on each slice. Remember to press Return when
% you have finished with that slice. You may skip ROI shapes by clicking
% outside the image field within the figure popup the same number of times
% as shapes you would like to skip (e.g., twice if you have
% drawn only one, or three if you draw no ROI shapes for that slice),
% and then press Return. 
%
% 4. Profit
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Housekeeping
    close all; 
    clear all;  
    clc; 
    
    %% Hard-coded inputs
    num_slices = 32; % Number of slices expected in volume 
    roi_shape_cutoff = 15; % Minimum number of voxels allowed for a contiguous ROI shape 
    roi_num_per_slice = 6; 

    %% Prepare and check input and output files and paths 
    % Prompt user for Nifti volume to mask 
    [file_to_load, file_to_load_dir] = uigetfile('*.nii'); 
    file_to_load_fullpath = fullfile(file_to_load_dir, file_to_load); 
    
    % Make sure selected file is in Nifti format
    [path, name, ext] = fileparts(file_to_load); 
    if strmatch(ext,'.nii') || strmatch(ext,'.nii.gz')
        disp("File confimed to be Nifti!") 
    end
    
    % Create output directory based on input file
    output_dir = sprintf('%sBrainMask', file_to_load_dir); 
    mkdir(output_dir); 
    
    % Prepare output file name and path based on input file
    output_file = sprintf('%s_VentricleMask_v2.nii', name)
    output_fullpath = fullfile(output_dir, output_file); 

    %% Load prechecked input file 
    % Load selected volume 
    vol_to_mask = niftiread(file_to_load_fullpath);
    vol_info = niftiinfo(file_to_load_fullpath); 

    vol_to_mask_rot = imrotate3(vol_to_mask, 270, [1 0 0], "nearest", "loose");
    vol_to_mask_rot = imrotate3(vol_to_mask_rot, 180, [0 0 1], "nearest", "loose"); 


    %% Prepare manually drawn ROIs
    % Initialize mask matrix based on input file
    mask_rot = zeros(size(vol_to_mask_rot)); 
    mask_rot_bin = mask_rot; 
    mask_rot_bin_clean = mask_rot; 
    
    %% Draw ROIs slice by slice and convert to 3D binary mask 
    % Walk through input volume slice-by-slice to draw slicewise mask
    for ii = 1:num_slices
        f1 = figure(); 
        colormap(gray)
        imagesc(squeeze(vol_to_mask_rot(ii,:,:))); 
        axis image; 
        set(gcf, 'Position', get(0, 'Screensize'));
            
        % Draw up to three ROIs per slice: Press Return when finished
        % adjusting
        for kk = 1:roi_num_per_slice
            outline(kk) = drawfreehand;
        end
        input(sprintf('Slice %d of %d: Press enter when done adjusting up to three ROIs', ii, num_slices));
    
        % Merge separate ROIs on same slice
        for kk = 1:roi_num_per_slice
            if kk == 1
                mask_outline_merged = createMask(outline(kk)); 
            else
                mask_outline_merged = mask_outline_merged + createMask(outline(kk)); 
            end
        end
    
        % Convert drawn ROIs to binary mask
        mask_bin_rot(ii, :, :) = imbinarize(mask_outline_merged); 
    
        % Clean up ROIs smaller than roi_shape_cutoff voxels 
        if any(any(mask_bin_rot(ii, :, :)))
            mask_bin_rot_clean(ii, :, :) = bwareafilt(logical(squeeze(mask_bin_rot(ii, :, :))),[roi_shape_cutoff Inf]); 
        else
            mask_bin_rot_clean(ii, :, :) = logical(squeeze(mask_bin_rot(ii, :, :))); 
        end
    
        close(f1); 
    
    end
    
    %% Rotate mask back to original orientation  
    mask_bin_clean = imrotate3(mask_bin_rot_clean, 180, [0 0 1], "nearest", "loose");
    mask_bin_clean = imrotate3(mask_bin_clean, 90, [1 0 0], "nearest", "loose");

    %% View final mask and write result to Nifti file for further processing 
    volumeViewer(vol_to_mask, int16(mask_bin_clean)); 
    niftiwrite(int16(mask_bin_clean),output_fullpath,vol_info); 

%end


