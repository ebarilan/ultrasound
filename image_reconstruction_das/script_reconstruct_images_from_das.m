function script_reconstruct_images_from_das(acquisition_type , phantom_type , data_type, IsEmbedded, numAngels)
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

% Input:
% acquisition_type:  1 = simulation || 2 = experiments
% phantom_typ:       1 = resolution & distorsion || 2 = contrast & speckle quality
% data_type:         1 = 'iq' || 2 = 'rf' || 3 = 'interp2FullFreqMap' || 4
% = 'interp2BBMapFreq' || 5 = 'idft' || 6 = 'nufft' || 7 = 'SPURS'


% clear all;
% close all;
home; % clc;
addpath(genpath('../src'));

if nargin == 3
    numAngels = 1;
    IsEmbedded = 1;
elseif nargin == 4
    numAngels = 1;
end

%-- Parameters
% acquisition_type = 1;       %-- 1 = simulation || 2 = experiments
% phantom_type = 1;           %-- 1 = resolution & distorsion || 2 = contrast & speckle quality
% data_type = 4;              %-- 1 = IQ || 2 = RF || 3 = IQInterpFFTfullFreqMap || 4
% = IQInterpFFTBaseband || 5 = IQInterpIDFTBaseband || 6 = nufft


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
    case {6,8}
        dataSave = 'nufft';
        
        run ../../../irt/setup.m
    case 7
        dataSave = 'SPURS';
        
        setenv('SPURS_DIR', '../../../SPURS_DEMO');
        SPURS_work_dir = getenv('SPURS_DIR');
        addpath(genpath(SPURS_work_dir));
    otherwise       %-- Do deal with bad values
        dataSave = 'iq';        
end


%-- Create path to load corresponding files
path_dataset = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_img_from_',dataSave,'numOfAngles',num2str(numAngels),'Is_embedad',num2str(IsEmbedded),'.hdf5'];
path_reconstruted_img_fig = ['../../reconstructed_image/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_img_from_',dataSave];

%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);


%-- Indices of plane waves to be used for each reconstruction
if numAngels == 1
    pw_indices{1} = 38;
else
    pw_indices{1} = round(linspace(1,dataset.firings,numAngels));
end
% pw_indices{1} = round(linspace(1,dataset.firings,11));
% pw_indices{1} = round(1:dataset.firings);               %-- dataset.firings corresponding to the total number of emitted steered plane waves

% IsEmbedded = 1;

%-- Reconstruct Bmode images for each pw_indices
disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
switch data_type    
    case 1
        image = das_iq(scan,dataset,pw_indices);
    case 2
        image = das_rf(scan,dataset,pw_indices);
    case {3,4,5,6,7}
        image = IQInterpFFT(scan,dataset,pw_indices, data_type, IsEmbedded);
    case 8 
        image = IQInterpFFT2(scan,dataset,pw_indices, data_type);
    otherwise       %-- Do deal with bad values
        image = das_iq(scan,dataset,pw_indices);       
end
disp('Reconstruction Done')
disp(['Result saved in "',path_reconstruted_img,'"'])


%-- Show the corresponding beamformed images
dynamic_range = 60;
[~,fig] = image.show(dynamic_range);
%-- Save results
savefig(fig,strcat(path_reconstruted_img_fig,'a',num2str(acquisition_type),'p', num2str(phantom_type),'t',num2str(data_type),'e',num2str(IsEmbedded),'a', num2str(numAngels),'_',datestr(now,'dd-mm-yy_HH-MM'),'.fig'));
saveas(fig,path_reconstruted_img_fig,'jpeg');
image.write_file(path_reconstruted_img);

if(0)
part = 0.55;
indNumPart = ceil(size(scan.z_matrix,1)*part);
[~,fig] = image.showPart(dynamic_range,[],indNumPart);
%-- Save results
savefig(fig,strcat(path_reconstruted_img_fig,datestr(now,'dd-mm-yy_HH-MM'),'_Part.fig'));
saveas(fig,strcat(path_reconstruted_img_fig,'_Part'),'jpeg');
end
