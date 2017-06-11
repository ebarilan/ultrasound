function [imageRec,imageFFT] = imageRecovery(acquisition_type , phantom_type , data_type, IsEmbedded, numAngels, spursConfig)

home; % clc;
addpath(genpath('../src'));

if nargin == 3
    numAngels = 1;
    IsEmbedded = 1;
elseif nargin == 4
    numAngels = 1;
end

%%
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

strSpurs = '';

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
        
        strSpurs = strcat('Sd',num2str(spursConfig.KernelFunctionDegree),...
                          'Si',num2str(spursConfig.Niterations),...
                          'So',num2str(spursConfig.OverGridFactor));
    case 9
        dataSave = 'intLinear1D';
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
if numAngels == 1
    pw_indices{1} = 38;
else
    pw_indices{1} = round(linspace(1,dataset.firings,numAngels));
end
% pw_indices{1} = round(linspace(1,dataset.firings,11));
% pw_indices{1} = round(1:dataset.firings);               %-- dataset.firings corresponding to the total number of emitted steered plane waves

% IsEmbedded = 1;
%%

%-- Reconstruct Bmode images for each pw_indices
disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
switch data_type    
    case 1
        image = das_iq(scan,dataset,pw_indices);
    case 2
        image = das_rf(scan,dataset,pw_indices);
    case {3,4,5,6,7,8,9}
        [image,imageFFT] = IQInterpFFT(scan,dataset,pw_indices, data_type, IsEmbedded, spursConfig);
    otherwise       %-- Do deal with bad values
        image = das_iq(scan,dataset,pw_indices);       
end
disp('Reconstruction Done')
disp(['Result saved in "',path_reconstruted_img,'"'])


%-- Show the corresponding beamformed images
dynamic_range = 60;
[~,fig] = image.show(dynamic_range);
%-- Save results

imageRec = image;