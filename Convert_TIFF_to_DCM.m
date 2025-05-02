function Convert_TIFF_to_DCM()

%% Directory and file handling

% User chooses TIFF directory
directory_TIF = uigetdir();

% Create a new directory for converted DCMs by modifying TIFF filepath
directory_DCM = strcat(directory_TIF, '_DCM_imposed_metadata'); 

% Move to TIFF directory
cd(directory_TIF);

% Define DCM info template: Need a sample DCM as input 
dcm_hdr = dicominfo('xxx.dcm'); 

% Find all tiff files in folder 
tiffStruct  = dir('*.tif');
tiffFileCell = {tiffStruct.name}'; 
tiffFileN = length(tiffFileCell); 

% Make new folder for dcm outputs 
mkdir(directory_DCM); 

%% Generate series-unique UIDs

StudyInstanceUID_info_new = dicomuid; 
SeriesInstanceUID_info_new = dicomuid; 
SliceLocation_info_new = 0; % Start at zero
SliceThickness_info_new = 0.004; % In mm  

ImagePositionPatient_x_info_new = 2048; 
ImagePositionPatient_y_info_new = 3891; % Note that y and z are flipped in order to force MICE Toolkit to recognize slices 
ImagePositionPatient_z_info_new = 0; % Note that y and z are flipped in order to force MICE Toolkit to recognize slices 

ImageType_info_new = 'DERIVED\PRIMARY\OTHER';

for ii=1:tiffFileN

    %% Convert TIFF to DCM
    % Locate TIFF file to convert in image series
    tiff_to_convert = tiffFileCell{ii};
    tiff_in_memory_big_uncropped = imread(tiff_to_convert); 

    % Crop the TIFF file from the top 
    % tiff_in_memory_big = imcrop(tiff_in_memory_big_uncropped,[2 1300 2046 2373]);

    % Time to downsample this enormous image
    % tiff_in_memory = imresize(tiff_in_memory_big, 1/5, 'bilinear');
    tiff_in_memory = imresize(tiff_in_memory_big_uncropped, 1/2, 'bilinear');
    
    % Set up new DCM filename and path 
    dcm_to_write = strrep(tiff_to_convert, '.tif', '.dcm'); 
    dcm_to_write_path = strcat(directory_DCM, '\\', dcm_to_write); 
    
    % Write TIFF image in memory into a DCM file
    % dicomwrite(tiff_in_memory, dcm_to_write_path, dcm_hdr); 
    dicomwrite(tiff_in_memory, dcm_to_write_path); 

    %% Change Patient ID and Study Date metadata %% 

    % Load in old dicom info
    DCM_TIFF_info_old = dicominfo(dcm_to_write_path); 

    % Determine appropriate Patient ID attribute from directory name
    %PatientID_info_folder_backslash_location = find(directory_TIF == '\', 1, 'last');
    %PatientID_info_new = extractAfter(directory_TIF, 137);
    PatientID_info_new = 'Light_Sheet_Dataset';

    % Dummy study date (should be refined) 
    StudyDate_info_new = '20230320';

    % Create new Study ID based solely on StudyDate (may change)
    StudyID_info_new = StudyDate_info_new;

    % Create new modality and image orientation 
    Modality_info_new = 'OT';
    PatientSpeciesDescription_info_new = 'RODENT';
    AnatomicalOrientationType_info_new = 'QUADRUPED';
    PatientPosition_info_new = 'HFS'; 
    ImagePositionPatient_info_new = [ImagePositionPatient_x_info_new; ImagePositionPatient_y_info_new; ImagePositionPatient_z_info_new]; 
    ImageOrientationPatient_info_new = [1; 0; 0; 0; 1; 0]; 
    SpacingBetweenSlices_info_new = 0; 

    % Determine location of Study Instance UID attribute in DCM struct and update
    StudyInstanceUID_info = dicomfind(DCM_TIFF_info_old,"StudyInstanceUID");
    StudyInstanceUID_info.Value = StudyInstanceUID_info_new;
    DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_old,StudyInstanceUID_info);

    % Determine location of Series Instance UID attribute in DCM struct and update
    SeriesInstanceUID_info = dicomfind(DCM_TIFF_info_old,"SeriesInstanceUID");
    SeriesInstanceUID_info.Value = SeriesInstanceUID_info_new;
    DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,SeriesInstanceUID_info);

    % Determine location of Patient ID attribute in DCM struct and update
    PatientID_info = dicomfind(DCM_TIFF_info_old,"PatientID");
    PatientID_info.Value = PatientID_info_new;
    DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,PatientID_info);

    % Determine location of Study Date attribute in DCM struct and update
    StudyDate_info = dicomfind(DCM_TIFF_info_old,"StudyDate");
    StudyDate_info.Value = StudyDate_info_new;
    DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,StudyDate_info);

    % % Determine location of Study ID attribute in DCM struct and update
    StudyID_info = dicomfind(DCM_TIFF_info_old,"StudyID");
    StudyID_info.Value = StudyID_info_new;
    DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,StudyID_info);

    %% Add series and instance information so individual files are viewed as a stack %%

    % Define new series and instance information based on data available 
    SeriesNumber_info_new = strrep(num2str(PatientID_info_new - 0), ' ', '');
    AcquisitionNumber_info_new = ii; 
    InstanceNumber_info_new = ii;
    ImagesInAcquisition_info_new = tiffFileN;

    % % Determine location of Series Number attribute in DCM struct and update
    SeriesNumber_info = dicomfind(DCM_TIFF_info_old,"SeriesNumber");
    SeriesNumber_info.Value = str2num(SeriesNumber_info_new(end-4:end));
    DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,SeriesNumber_info);
 
    % Change modality from "Other" to "General Microscopy"
    Modality_info = dicomfind(DCM_TIFF_info_old,"Modality");
    Modality_info.Value = Modality_info_new;
    DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,Modality_info);

    % Determine location of Image Type attribute in DCM struct and update
    try
        ImageType_info = dicomfind(DCM_TIFF_info_old,"ImageType");
        ImageType_info.Value = ImageType_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,ImageType_info);
    catch
        DCM_TIFF_info_new.ImageType = ImageType_info_new;
    end

    % Determine location of Patient Species Description attribute in DCM struct and update
    try
        PatientSpeciesDescription_info = dicomfind(DCM_TIFF_info_old,"PatientSpeciesDescription");
        PatientSpeciesDescription_info.Value = PatientSpeciesDescription_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,PatientSpeciesDescription_info);
    catch
        DCM_TIFF_info_new.PatientSpeciesDescription = PatientSpeciesDescription_info_new;
    end

    % Determine location of Anatomical Orientation Type attribute in DCM struct and update
    try
        AnatomicalOrientationType_info = dicomfind(DCM_TIFF_info_old,"AnatomicalOrientationType");
        AnatomicalOrientationType_info.Value = AnatomicalOrientationType_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,AnatomicalOrientationType_info);
    catch
        DCM_TIFF_info_new.AnatomicalOrientationType = AnatomicalOrientationType_info_new;
    end

    % Determine location of Patient Position attribute in DCM struct and update
    try
        PatientPosition_info = dicomfind(DCM_TIFF_info_old,"PatientPosition");
        PatientPosition_info.Value = PatientPosition_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,PatientPosition_info);
    catch
        DCM_TIFF_info_new.PatientPosition = PatientPosition_info_new;
    end

    % Determine location of Image Position Patient attribute in DCM struct and update
    try
        ImagePositionPatient_info = dicomfind(DCM_TIFF_info_old,"ImagePositionPatient");
        ImagePositionPatient_info.Value = ImagePositionPatient_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,ImagePositionPatient_info);
    catch
        DCM_TIFF_info_new.ImagePositionPatient = ImagePositionPatient_info_new;
    end

    % Determine location of Image Orientation Patient attribute in DCM struct and update
    try
        ImageOrientationPatient_info = dicomfind(DCM_TIFF_info_old,"ImageOrientationPatient");
        ImageOrientationPatient_info.Value = ImageOrientationPatient_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,ImageOrientationPatient_info);
    catch
        DCM_TIFF_info_new.ImageOrientationPatient = ImageOrientationPatient_info_new;
    end

    % Determine location of Acquisition Number attribute in DCM struct and update
    try
        AcquisitionNumber_info = dicomfind(DCM_TIFF_info_old,"AcquisitionNumber");
        AcquisitionNumber_info.Value = AcquisitionNumber_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,AcquisitionNumber_info);
    catch
        DCM_TIFF_info_new.AcquisitionNumber = AcquisitionNumber_info_new;
    end

    %Determine location of Instance Number attribute in DCM struct and update
    try
        InstanceNumber_info = dicomfind(DCM_TIFF_info_old,"InstanceNumber");
        InstanceNumber_info.Value = InstanceNumber_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,InstanceNumber_info);
    catch
        DCM_TIFF_info_new.InstanceNumber = InstanceNumber_info_new;
    end

    % Determine location of Images in Acquisition attribute in DCM struct and update
    try
        ImagesInAcquisition_info = dicomfind(DCM_TIFF_info_old,"ImagesInAcquisition");
        ImagesInAcquisition_info.Value = ImagesInAcquisition_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,ImagesInAcquisition_info);
    catch
        DCM_TIFF_info_new.ImagesInAcquisition = ImagesInAcquisition_info_new;
    end

    %Determine location of Slice Location attribute in DCM struct and update
    try
        SliceLocation_info = dicomfind(DCM_TIFF_info_old,"SliceLocation");
        SliceLocation_info.Value = SliceLocation_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,SliceLocation_info);
    catch
        DCM_TIFF_info_new.SliceLocation = SliceLocation_info_new;
    end

    %Determine location of Slice Thickness attribute in DCM struct and update
    try
        SliceThickness_info = dicomfind(DCM_TIFF_info_old,"SliceThickness");
        SliceThickness_info.Value = SliceThickness_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,SliceThickness_info);
    catch
        DCM_TIFF_info_new.SliceThickness = SliceThickness_info_new;
    end

    %Determine location of Spacing Between Slices attribute in DCM struct and update
    try
        SpacingBetweenSlices_info = dicomfind(DCM_TIFF_info_old,"SpacingBetweenSlices");
        SpacingBetweenSlices_info.Value = SpacingBetweenSlices_info_new;
        DCM_TIFF_info_new = dicomupdate(DCM_TIFF_info_new,SpacingBetweenSlices_info);
    catch
        DCM_TIFF_info_new.SpacingBetweenSlices = SpacingBetweenSlices_info_new;
    end

    %% Write new DCM including Patient ID, Study Date, and series information %%

    dicomwrite(tiff_in_memory, dcm_to_write_path, DCM_TIFF_info_new, "CreateMode", "copy");  
    
    % Report progress 
    fprintf('Slice %d TIFF converted to DCM!\n', ii)

    % Update slice location
    SliceLocation_info_new = SliceLocation_info_new + SliceThickness_info_new;
    ImagePositionPatient_z_info_new = ImagePositionPatient_z_info_new + 1;

end

end