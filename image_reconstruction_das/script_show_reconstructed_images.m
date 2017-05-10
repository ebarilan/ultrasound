%-- Script to be used as an example to display the reconstructed images

%-- After choosing the specific configuration through acquisition_type, 
%-- phantom_type and data_type parameters, this script allows displaying the 
%-- reconstructed iamges (by the default, the proposed implemented technique 
%-- corresponds to the standard Delay And Sum (DAS) technique with apodization in reception

%-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
%--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)

%-- $Date: 2016/03/01 $  


clear all;
close all;
clc;
addpath(genpath('../src'));


%-- Parameters
acquisition_type = 1;       %-- 1 = simulation || 2 = experiments
phantom_type = 1;           %-- 1 = resolution || 2 = contrast
data_type = 1;              %-- 1 = IQ || 2 = RF


%-- Parsing parameter choices
switch acquisition_type    
    case 1
        acquisition = 'simulation';
        acqui = 'simu';
    case 2
        acquisition = 'experiments';
        acqui = 'expe';
    otherwise       %-- Do deal with bad values
        acquisition = 'simulation';
        acqui = 'simu';        
end
switch phantom_type    
    case 1
        phantom = 'resolution_distorsion';
    case 2
        phantom = 'contrast_speckle';
    otherwise       %-- Do deal with bad values
        phantom = 'resolution';
end
switch data_type    
    case 1
        data = 'iq';
    case 2
        data = 'rf';
    otherwise       %-- Do deal with bad values
        data = 'iq';        
end


%-- Create path to load corresponding files
path_dataset = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_img_from_',data,'.hdf5'];


%-- If corresponding image file exists then display it
if exist(path_reconstruted_img,'file')

    %-- Read corresponding files
    disp(['Show reconstructed image from ',acquisition,' for ',phantom,' using ',data,' dataset'])
    dataset = us_dataset();
    dataset.read_file(path_dataset);
    scan = linear_scan();
    scan.read_file(path_scan);
    image = us_image();
    image.read_file(path_reconstruted_img);

    %-- Show the beamformed images
    dynamic_range = 60;
    image.show(dynamic_range);

else 
    disp(['File: "',path_reconstruted_img,'" does not exists'])
end

