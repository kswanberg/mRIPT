% Whole-brain MRI time series analysis 
% 
% Loads two folders of time-series magnetic resonance
% image DICOMs, converts them to NII, registers each volume to the first
% time point, applies a skull-stripping mask expected to be NIFTI format,
% and combines all registered, skull-stripped volumes into a single animated GIF
%
% Dependencies: Requires Xiang Ruili's dicm2nii
% (https://github.com/xiangruili/dicm2nii) to be in MATLAB path
% 
% Inputs: 
% 
%  1. Time-series MRI volumes split into two folders (can adjust number of
%  expected volumes in input list) 
% 
%  2. Skull-stripped NIFTI mask corresponding to first time-series volume
%  (i.e. time point 1)
% 
% Outputs: 
% 
% 1. Unprocessed NIFTI files for each volume
% 
% 2. Registered (to first time point) NIFTI files for each volume
% 
% 3. Registered and skull-stripped NIFTI files for each volume 
% 
% 4. Animated GIF of all registered and skull-stripped time points. Note
% that the animated GIF labeling starting around line 270 is specific to a
% DCE-MRI series begun 40 minutes after contrast injection and can be
% manually adjusted for other applications. 
% 
% Author: Kelley Swanberg (Lund University, 2024) 
% 
% 


%% Housekeeping
close all; 
clear all; 
clc; 

%% Inputs 
% Filepaths 
exp_id = 'IP_01'; 
folder_1 = '191121_IL_MR94_glymphatic_system_20240308_Efflux___E3_P1'; % 15 minutes
folder_2 = '191121_IL_MR94_glymphatic_system_20240308_Efflux___E11_P1'; % 50 minutes
be_mask = 'Efflux_1mcLpmin_38_t1\wholemask.nii'; % BE mask 

% Other variables
num_slices = 100; % number of slices per time point 
folder_1_num_time_points = 6; 
folder_2_num_time_points = 20;
folder_1_num_dcms_expected = folder_1_num_time_points*num_slices; 
folder_2_num_dcms_expected = folder_2_num_time_points*num_slices; 
volume_acquisition_duration = 2.5 % minutes
total_time_points = folder_1_num_time_points + folder_2_num_time_points; 
nii_master_folder = sprintf('%s_nii_time_series', exp_id)
nii_master_folder_reg = sprintf('%s_nii_time_series_reg', exp_id)
nii_master_folder_reg_ss = sprintf('%s_nii_time_series_reg_ss', exp_id)

%% Algorithm 

%% Sort DCMs into individual folders
% Clean up folder 1 path info
folder_1_pathinfo = dir(fullfile(folder_1, '*.dcm'));

% Clean up folder 2 path info
folder_2_pathinfo = dir(fullfile(folder_2, '*.dcm'));

% Check that folder 1 has expected number of files and warn user if not the
% case
folder_1_dcmarraysize = size(folder_1_pathinfo); 
folder_1_num_dcms = folder_1_dcmarraysize(1);

if folder_1_num_dcms == folder_1_num_dcms_expected
    disp("Folder 1 file number check passed!") 
end

% Check that folder 2 has expected number of files and warn user if not the
% case
folder_2_dcmarraysize = size(folder_2_pathinfo); 
folder_2_num_dcms = folder_2_dcmarraysize(1);

if folder_2_num_dcms == folder_2_num_dcms_expected
    disp("Folder 2 file number check passed!")
end

% Sort dicoms into individual folders for time series analysis 
% Start with folder 1
for ii = 1:folder_1_num_time_points
    folder_name = sprintf('%s_%d', exp_id, ii); 
    mkdir(folder_name);
    first_slice = (ii - 1)*num_slices + 1; 

    for kk = 1:num_slices
        slice_to_copy = first_slice + (kk - 1); 
        file_name = folder_1_pathinfo(slice_to_copy).name; 
        file_source = sprintf('%s/%s', folder_1, file_name); 
        file_destination = folder_name; 
        copyfile(file_source, file_destination );
    end
end

% Continue with folder 2
folder_2_max_time_points = folder_1_num_time_points + folder_1_num_time_points; 

for ii = 1:folder_2_num_time_points
    folder_2_tp = folder_1_num_time_points + ii; 
    folder_name = sprintf('%s_%d', exp_id, folder_2_tp); 
    mkdir(folder_name);
    first_slice = (ii - 1)*num_slices + 1; 

    for kk = 1:num_slices
        slice_to_copy = first_slice + (kk - 1); 
        file_name = folder_2_pathinfo(slice_to_copy).name; 
        file_source = sprintf('%s/%s', folder_2, file_name); 
        file_destination = folder_name; 
        copyfile(file_source, file_destination);
    end
end

mkdir(nii_master_folder); 
mkdir(nii_master_folder_reg); 
mkdir(nii_master_folder_reg_ss); 

%% Convert DCMs in individual folders into .nii
for ii = 1:total_time_points
    folder_name = sprintf('%s_%d', exp_id, ii); 
    folder_name_nii = sprintf('%s_nii', folder_name); 
    dicm2nii(folder_name, folder_name_nii, 0); 

    if ii == 1
        baseline_folder_name_nii = folder_name_nii; 
    end

    % Rename all files in folder to folder name
    nii_folder_pathinfo = dir(folder_name_nii);
    nii_folder_pathinfo_clean = nii_folder_pathinfo(3:end); 

    nii_folder_arraysize = size(nii_folder_pathinfo_clean); 
    nii_folder_num_files = nii_folder_arraysize(1);

    for kk = 1:nii_folder_num_files 
        old_file = sprintf('%s/%s', folder_name_nii, nii_folder_pathinfo_clean(kk).name); 
        [filepath,name,ext] = fileparts(old_file); 
        new_file = sprintf('%s/%s%s', filepath, folder_name_nii, ext); 
        copyfile(old_file, new_file); 
        delete(old_file); 

        new_file_master_folder = sprintf('%s/%s%s', nii_master_folder, folder_name_nii, ext); 
        copyfile(new_file, new_file_master_folder); 
    end
end

%% Register all .nii to baseline 
baseline_volume_file = sprintf('%s/%s.nii', baseline_folder_name_nii, baseline_folder_name_nii); 
baseline_volume = niftiread(baseline_volume_file);

% Set up registration defaults 
[optimizer,metric] = imregconfig("multimodal"); % Because DCE-MRI designs can dramatically affect image contrast over time
optimizer.InitialRadius = 0.004;

for ii = 2:total_time_points
    folder_name = sprintf('%s_%d', exp_id, ii); 
    folder_name_nii = sprintf('%s_nii', folder_name); 

    volume_file = sprintf('%s/%s.nii', folder_name_nii, folder_name_nii); 
    volume_file_volume_info = niftiinfo(volume_file);
    volume_file_volume = niftiread(volume_file);

    % Register volume to baseline while ensuring consistent matrix size from input to output 
    reg_tform = imregtform(volume_file_volume,baseline_volume,"rigid",optimizer,metric);
    reg_tform_R = affineOutputView(size(volume_file_volume),reg_tform,"BoundsStyle","SameAsInput");
    [volume_reg, volume_file_reg_ref] = imwarp(volume_file_volume, reg_tform, "OutputView", reg_tform_R); 

    % volume_file_registered = imregister(volume_file_volume, baseline_volume, "rigid",optimizer,metric);
    %Save to same folder
    volume_file_reg = strrep(volume_file, '.nii', '_reg.nii'); 
    niftiwrite(volume_reg,volume_file_reg,volume_file_volume_info); 

    % Also save to master folder
    volume_file_reg_master = sprintf('%s/%s_reg.nii', nii_master_folder_reg, folder_name_nii);
    niftiwrite(volume_reg,volume_file_reg_master,volume_file_volume_info)
end

%% Skull-strip all registered .nii using existing BE mask
mask_file_volume = niftiread(be_mask);
mask_file_volume_flipped = flipud(mask_file_volume); % Because this is upside-down relative to brains

for ii = 1:total_time_points
    folder_name = sprintf('%s_%d', exp_id, ii); 
    folder_name_nii = sprintf('%s_nii', folder_name); 
    
    if ii == 1
        volume_file_reg = sprintf('%s/%s.nii', folder_name_nii, folder_name_nii); 
    else
        volume_file_reg = sprintf('%s/%s_reg.nii', folder_name_nii, folder_name_nii); 
    end

    volume_file_volume_info_reg= niftiinfo(volume_file_reg);
    volume_file_volume_reg = niftiread(volume_file_reg);

    % Apply mask 
    volume_reg_ss = volume_file_volume_reg .* int16(imbinarize(mask_file_volume_flipped)); 

    %Save to same folder
    volume_file_reg_ss = strrep(volume_file_reg, '.nii', '_ss.nii'); 
    niftiwrite(volume_reg_ss,volume_file_reg_ss,volume_file_volume_info_reg); 

    % Also save to master folder
    volume_file_reg_ss_master = sprintf('%s/%s_ss.nii', nii_master_folder_reg_ss, folder_name_nii);
    niftiwrite(volume_reg_ss,volume_file_reg_ss_master,volume_file_volume_info_reg)
end

%% View time series and export animated GIF
viewer1 = viewer3d(BackgroundColor="black", BackgroundGradient = "off"); 
viewer1.CameraPosition = [118.2719 150.9215 10.1190]; 
viewer1.CameraTarget = [45.8359 69.6273 87.8687]; 
viewer1.CameraUpVector = [0.5579 -0.4233 0.7139];
%viewer1.CameraZoom = 1.7969; 
viewer1.CameraZoom = 1.3; % For Nte
viewer1.OrientationAxes = 'off'; 

viewer2 = viewer3d(BackgroundColor="black", BackgroundGradient = "off"); 
viewer2.CameraPosition = [179.0412 76.1435 73.1783]; 
viewer2.CameraTarget = [45.9829 79.6222 86.7482];
viewer2.CameraUpVector = [-0.2736 0.0122 0.9618];
%viewer2.CameraZoom = 1.4324;
viewer2.CameraZoom = 1; % For Nte
viewer2.OrientationAxes = 'off'; 

viewer3 = viewer3d(BackgroundColor="black", BackgroundGradient = "off"); 
viewer3.CameraPosition = [54.9178 208.3474 49.0443];
viewer3.CameraTarget = [47.5936 78.4419 80.2159];
viewer3.CameraUpVector = [0.9052 -0.4235 0.0357];
%viewer3.CameraZoom = 1.6800;
viewer3.CameraZoom = 1.3; % For Nte
viewer3.OrientationAxes = 'off'; 

viewer4 = viewer3d(BackgroundColor="black", BackgroundGradient = "off"); 
viewer4.CameraPosition = [5.9477 81.8555 226.1828]; 
viewer4.CameraTarget = [38.3513 79.1536 96.4005]; 
viewer4.CameraUpVector = [0.9904 -0.0384 -0.1334];
%viewer4.CameraZoom = 2.5453;
viewer4.CameraZoom = 2.3; % For Nte
viewer4.OrientationAxes = 'off'; 

filename_gif = sprintf("%s_reg_ss_anim.gif", exp_id); % Specify the output file name

for ii = 1:total_time_points
    filename = sprintf('%s/%s_%d_nii_ss.nii',nii_master_folder_reg_ss, exp_id, ii);
    v = niftiread(filename); 
        if ii == 1
            baseline = v; 
        end
    vol1 = volshow(v, RenderingStyle="MaximumIntensityProjection",  Parent=viewer1);
    pause(3); 
    vol2 = volshow(v, RenderingStyle="MaximumIntensityProjection",  Parent=viewer2);
    pause(3); 
    vol3 = volshow(v, RenderingStyle="MaximumIntensityProjection",  Parent=viewer3);
    pause(3); 
    vol4 = volshow(v, RenderingStyle="MaximumIntensityProjection",  Parent=viewer4);
    pause(3); 

    % Adapted from https://se.mathworks.com/matlabcentral/answers/
    % 1961729-how-to-save-volshow-data-config-as-fig-or-img-volumeviewer
    obj1 = ancestor(vol1,'figure','toplevel');
    F1 = getframe(obj1);
    
    obj2 = ancestor(vol2,'figure','toplevel');
    F2 = getframe(obj2);

    obj3 = ancestor(vol3,'figure','toplevel');
    F3 = getframe(obj3);

    obj4 = ancestor(vol4,'figure','toplevel');
    F4 = getframe(obj4);

    figure();
    
    % Note that this time labeling is specific to a DCE-MRI series begun 40 minutes after contrast injection
    minutes_post_a = round(40 + (ii-1)*volume_acquisition_duration);  
    minutes_post_b = round(40 + (ii)*volume_acquisition_duration);
    row1 = horzcat(F1.cdata, F2.cdata); 
    row2 = horzcat(F3.cdata, F4.cdata); 
    whole = vertcat(row1, row2); 
    fullgrid = imagesc(whole);
    titlecap = sprintf('%.0f-%.0f minutes after injection', minutes_post_a, minutes_post_b); 
    h = title(titlecap); 
    set(gcf, 'color', 'black');
    set(gcf, 'InvertHardCopy', 'off');
    set(h, 'Color','white'); 
    grid off; 
    axis off; 

    filename_png = strrep(filename, '.nii', '.png'); 
    saveas(fullgrid, filename_png, 'png'); 

    % Animated GIF
    [A,map] = rgb2ind(whole,256);
    if ii == 1
       imwrite(A,map,filename_gif,"gif","LoopCount",Inf,"DelayTime",0.18);
    else
       imwrite(A,map,filename_gif,"gif","WriteMode","append","DelayTime",0.18);
    end
    close all; 
end



