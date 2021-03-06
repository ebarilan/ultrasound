function script_reconstruct_images_from_das_TEST_Tanya(acquisition_type , phantom_type , data_type)
%-- Script to be used as an example to manipulate the provided dataset

%-- After choosing the specific configuration through acquisition_type,
%-- phantom_type and data_type parameters, this script allows reconstructing
%-- images for evaluation for different choices of steered plane waves involved
%-- in the compounding scheme (specified by the pw_indices parameter)

%-- The implemented method (to be used as example) corresponds to the standard Delay
%-- And Sum (DAS) technique with apodization in reception

%-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
%--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)

%-- $Date: 2016/03/01 $


% clear all;
% close all;
home; % clc;
addpath(genpath('../src'));


%-- Parameters
acquisition_type = 1;       %-- 1 = simulation || 2 = experiments
phantom_type = 1;           %-- 1 = resolution & distorsion || 2 = contrast & speckle quality
data_type = 4;              %-- 1 = IQ || 2 = RF || 3 = IQInterpFFTfullFreqMap || 4
% = IQInterpFFTBaseband || 5 = IQInterpIDFTBaseband || 6 = nufft || 7 = padded data X

pulseFlag = 0; %-- 0 = pulse shape is not considered || 1 = pulse shape is computed

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

switch data_type
    case 1
        dataSave = 'iq';
    case 2
        dataSave = 'rf';
    case 3
        dataSave = 'interp2FullFreqMap';
    case 4
        dataSave = 'interp2BBMapFreq';
    case 5
        dataSave = 'idft';
    case 6
        dataSave = 'nufft';
    case 7
        dataSave = 'paddedX';
    case 8
        dataSave = 'paddedZ';
    otherwise       %-- Do deal with bad values
        
        dataSave = 'iq';
end


%-- Create path to load corresponding files
path_dataset = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_img_from_',dataSave,'.hdf5'];
path_reconstruted_img_fig = ['../../reconstructed_image/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_img_from_',dataSave];


%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);


%-- Indices of plane waves to be used for each reconstruction
pw_indices{1} = 38;
% pw_indices{2} = round(linspace(1,dataset.firings,3));
% pw_indices{3} = round(linspace(1,dataset.firings,11));
% pw_indices{4} = round(1:dataset.firings);               %-- dataset.firings corresponding to the total number of emitted steered plane waves


%-- Reconstruct Bmode images for each pw_indices
disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
switch data_type
    case 1
        image = das_iq(scan,dataset,pw_indices);
    case 2
        image = das_rf(scan,dataset,pw_indices);
    case {3,4,5,6,7,8}
        image = IQInterpFFT(scan,dataset,pw_indices, data_type,pulseFlag);
    otherwise       %-- Do deal with bad values
        image = das_iq(scan,dataset,pw_indices);
end
disp('Reconstruction Done')
disp(['Result saved in "',path_reconstruted_img,'"'])


%-- Show the corresponding beamformed images
dynamic_range = 60;
[~,fig] = image.show(dynamic_range);

%-- Save results
savefig(fig,sprintf(path_reconstruted_img_fig,'.fig'));
saveas(fig,path_reconstruted_img_fig,'jpeg');
image.write_file(path_reconstruted_img);

