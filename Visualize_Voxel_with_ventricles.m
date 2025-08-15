% Show 1H-MRS ROI on DCM

% Read stack of DCMs in folder
selpath = uigetdir; % Directory 
pathinfo = dir(selpath);
% Clean up path info
pathinfo_clean = pathinfo(3:end, :);

% Initialize stack of DCMs 
dcmarraysize = size(pathinfo_clean); 
num_dcms = dcmarraysize(1);

% Load stack of DCMs
for ii=1:num_dcms
dcmpath = sprintf('%s\\%s', pathinfo_clean(ii).folder, pathinfo_clean(ii).name);

% Store minimum and maximum pixel values 
dcm_minpix = repmat(0, [1, num_dcms]);
dcm_maxpix = repmat(0, [1, num_dcms]);

% If we are loading the first DCM we have to initialize our array 
if ii==1
    init_info = dicominfo(dcmpath); 
    init_dcm = dicomread(dcmpath); 
    dcm_pixel_samples = init_info.SamplesPerPixel;

    dcm_dim = size(init_dcm); 
    dcm_x_dim = dcm_dim(1); 
    dcm_y_dim = dcm_dim(2);
    dcm_stack = repmat(int16(0), [dcm_x_dim dcm_y_dim dcm_pixel_samples num_dcms]);
    dcm_stack_info = repmat(init_info, 1, num_dcms); 
end
    dcm_slice_info = dicominfo(dcmpath); 
    dcm_slice = dicomread(dcmpath); 
    dcm_stack_info(1, ii) = dcm_slice_info; 
    
    for jj = 1:dcm_pixel_samples
        dcm_stack(:, :, jj, ii) = dcm_slice; 
        dcm_minpix(ii) = dcm_slice_info.SmallestImagePixelValue; 
        dcm_maxpix(ii) = dcm_slice_info.LargestImagePixelValue; 
    end
end

% Rescale image to start at 0
% (https://se.mathworks.com/company/newsletters/articles/accessing-data-in-dicom-files.html)
b = min(dcm_minpix);
m = 2^16/(max(dcm_maxpix) - b);
dcm_stack_rescaled = imlincomb(m, dcm_stack, -(m * b), 'uint16');

load mri
dcm_stack_rescaled_sq = squeeze(dcm_stack_rescaled); 
figure

dcm_stack_rescaled_sq_idx = cell(dcm_x_dim, dcm_y_dim, num_dcms); 
dcm_stack_rescaled_sq_coord = dcm_stack_rescaled_sq_idx; 

% Origin from position information 
x_pos = init_info.ImagePositionPatient(1); % -8.7948
y_pos = init_info.ImagePositionPatient(2); % -5.6764
z_pos = init_info.ImagePositionPatient(3); % 12.5060

vx_spacing_z = init_info.PixelSpacing(1); % 0.1mm
vx_spacing_x = init_info.PixelSpacing(2); % 0.1mm
vx_spacing_y = init_info.SliceThickness; % 0.1mm

% But order of slices is z, x, y
% Create coordinate system matrix 
coord_matrix = [-vx_spacing_z 0 0 z_pos; 0 vx_spacing_x 0 x_pos; 0 0 vx_spacing_y y_pos; 0 0 0 1]; % Switch y and z here 

for ii = 1:dcm_x_dim
    for jj = 1:dcm_y_dim
        for kk = 1:num_dcms
            a =  [ii-1, jj-1, kk-1, 1]; 
            dcm_stack_rescaled_sq_idx{ii, jj, kk} = a'; 
            dcm_stack_rescaled_sq_coord{ii, jj, kk} = coord_matrix*a'; 
        end
    end
end

% Create voxel mask. Voxel provided in RPS coordinates 
voxel_center_x = 1.17318; % Negate from method file  
voxel_center_y = -1.4; % Negate from method file 
voxel_center_z = 3.93855;

voxel_dim_x = 2; 
voxel_dim_y = 2; 
voxel_dim_z = 2; 
 
% RPS; appears to be the center of the voxel 
voxel_max_x = voxel_center_x + 0.5*voxel_dim_x - vx_spacing_x; 
voxel_min_x = voxel_center_x - 0.5*voxel_dim_x; 
%voxel_max_x = voxel_center_x + 0.5*voxel_dim_x; 
%voxel_min_x = voxel_center_x - 0.5*voxel_dim_x +  vx_spacing_x; 
%voxel_min_y = voxel_center_y - 0.5*voxel_dim_y; 
voxel_max_y = voxel_center_y + 0.5*voxel_dim_y; 
%voxel_max_y = voxel_center_y + 0.5*voxel_dim_y - vx_spacing_y; 
%voxel_min_y = voxel_center_y - 0.5*voxel_dim_y; 
voxel_min_y = voxel_center_y - 0.5*voxel_dim_y + vx_spacing_y; 
%voxel_min_z = voxel_center_z - 1*voxel_dim_z; 
%voxel_max_z = voxel_center_z + 0*voxel_dim_z; 
voxel_min_z = voxel_center_z - 0.5*voxel_dim_z + vx_spacing_z; 
voxel_max_z = voxel_center_z + 0.5*voxel_dim_z; 
%voxel_min_z = voxel_center_z - 0.5*voxel_dim_z; 
%voxel_max_z = voxel_center_z + 0.5*voxel_dim_z - vx_spacing_z; 

% Remember that order of matrix is z, x, y
coord_x_vector_array = [dcm_stack_rescaled_sq_coord{1, :, 1}]; 
coord_x_vector = coord_x_vector_array(2, :); 
[d_x, ix_x] = min(abs(coord_x_vector - voxel_center_x));
[d_xmin, ix_xmin] = min(abs(coord_x_vector - voxel_min_x));
[d_xmax, ix_xmax] = min(abs(coord_x_vector - voxel_max_x));

coord_y_vector_array = [dcm_stack_rescaled_sq_coord{1, 1, :}]; 
coord_y_vector = coord_y_vector_array(3, :); 
[d_y, ix_y] = min(abs(coord_y_vector - voxel_center_y));
[d_ymin, ix_ymin] = min(abs(coord_y_vector - voxel_min_y));
[d_ymax, ix_ymax] = min(abs(coord_y_vector - voxel_max_y));

coord_z_vector_array = [dcm_stack_rescaled_sq_coord{:, 1, 1}]; 
coord_z_vector = coord_z_vector_array(1, :); 
[d_z, ix_z] = min(abs(coord_z_vector - voxel_center_z));
[d_zmin, ix_zmin] = min(abs(coord_z_vector - voxel_min_z));
[d_zmax, ix_zmax] = min(abs(coord_z_vector - voxel_max_z));

dcm_stack_voxel_mask = repmat(uint8(0), [dcm_x_dim dcm_y_dim num_dcms]);

% Create voxel mask. 
for ii = 1:dcm_x_dim % Actually zdim 
    for jj = 1:dcm_y_dim % Actually xdim 
        for kk = 1:num_dcms % Actually ydim

            px_bin = 0; 
            
            if ii <= ix_zmin && ii >= ix_zmax % had to flip this because z is backward
                if jj >= ix_xmin && jj <= ix_xmax
                    if kk >= ix_ymin && kk <= ix_ymax
                        px_bin = 1;  
                    end
                end
            end
            
            edgelord = 0; 
            if ii == ix_zmin || ii == ix_zmax % had to flip this because z is backward
    	        edgelord=edgelord+1; 
            end

            if jj == ix_xmin || jj == ix_xmax % had to flip this because z is backward
    	        edgelord=edgelord+1; 
            end

            if kk == ix_ymin || kk == ix_ymax % had to flip this because z is backward
    	        edgelord=edgelord+1; 
            end

            if edgelord >= 2 && px_bin == 1;  
                px_bin = 2; 
            end

            dcm_stack_voxel_mask(ii, jj, kk) = px_bin; 

        end
    end
end

% Initial conditions for volume viewer
alpha = [0 0 0.3 0.9];
color = [0 0 0; 200 140 75; 231 208 141; 255 255 255]./255;
intensity = [1000 2000 3000 4000];
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);

viewer = viewer3d(BackgroundColor="black",BackgroundGradient="off");

    %'VolumeRendering'
    %'MaximumIntensityProjection'
    %'MinimumIntensityProjection'
    %'Isosurface'
    %'SlicePlanes'
    %'GradientOpacity'
    %'CinematicRendering'
    %'LightScattering'

%image_and_mask = int16(dcm_stack_rescaled_sq) + dcm_stack_voxel_mask; 

% Non-skull-stripped
brainmask = niftiread('E:\Lactate_ISMRM_2025\Exp_2\Skull_Stripping\Data_to_Deidentify\191121_IL_MR94_glymphatic_system_20241109_iso_kx_fMRS_01_35587063\nii\brain_mask.nii.gz');
brainmask_rot1 = imrotate3(brainmask, 90, [0 -1 0], "nearest", "loose");
brainmask_rot2 = fliplr(uint16(imbinarize(imrotate3(brainmask_rot1, 90, [0 0 1], "nearest", "loose"))));
brainmask_flip = fliplr(flip(brainmask_rot2, 3)); 

masked_dcm = brainmask_rot2 .* dcm_stack_rescaled_sq; 

voxelmask_flip = fliplr(flip(dcm_stack_voxel_mask, 3)); 

ventmask = niftiread('E:\Lactate_ISMRM_2025\Exp_2\Skull_Stripping\Data_to_Deidentify\191121_IL_MR94_glymphatic_system_20241109_iso_kx_fMRS_01_35587063\nii\ventricle_mask.nii.gz');
ventmask_rot1 = imrotate3(ventmask, 90, [0 -1 0], "nearest", "loose");
ventmask_rot2 = fliplr(uint16(imbinarize(imrotate3(ventmask_rot1, 90, [0 0 1], "nearest", "loose"))));
ventmask_flip = fliplr(flip(ventmask_rot2, 3)); 

voxel_and_vent = voxelmask_flip + 3*uint8(ventmask_flip); 

%dcm_stack_rescaled_sq_flipped = flip(dcm_stack_rescaled_sq, 3);
%volumeViewer(fliplr(dcm_stack_rescaled_sq), fliplr(dcm_stack_voxel_mask)); % Because DCMs loaded here have a lr flip relative to those shown in PV
%volumeViewer(fliplr(flip(masked_dcm, 3)), fliplr(flip(dcm_stack_voxel_mask, 3)));
% volumeViewer(fliplr(flip(dcm_stack_rescaled_sq, 3)), fliplr(flip(dcm_stack_voxel_mask, 3)));

volumeViewer(fliplr(flip(masked_dcm, 3)), voxel_and_vent);

%volumeViewer(fliplr(flip(dcm_stack_rescaled_sq, 3)), ventmask_flip);


%volumeViewer(dcm_stack_rescaled_sq, dcm_stack_voxel_mask);
% 
% nii_stack_rescaled_sq_rot2_flipped = flip(nii_stack_rescaled_sq_rot2, 3); 

% 
% %Skull-stripped
% volumeViewer(nii_stack_rescaled_sq_rot2_flipped, dcm_stack_voxel_mask_flipped);

% Export voxel slicewise within image coordinate system as DICOM

% for ii=1:100
%     filename_dcm = sprintf("%s\\Voxel_slice%d.dcm", fileparts(dcmpath), ii); 
%     dicomwrite(int16(dcm_stack_voxel_mask(:, :, ii)),filename_dcm, dcm_stack_info(1, ii), 'CreateMode','Copy');
% end

% vol = volshow(dcm_stack_rescaled_sq, ...
%               Parent=viewer, ...
%               RenderingStyle='VolumeRendering', ...
%               Alphamap=alphamap, overlayData=dcm_stack_voxel_mask);