acquisition_type = 1;
ImageOpt = [1,2];

if(0)
%% 0.NUFFT
ProccesOpt = 6;
numAngels = [11];% numAngels = [1,3,5,7,9,11];
IsEmbedded = 1;%[0,1];
runScripts(numAngels, ProccesOpt, IsEmbedded, ImageOpt, acquisition_type, []);
end
%% 1.Spurs
ProccesOpt = 7;

spursConfig.Rho = 1e-3;
spursConfig.FilterInImageSpace = 1;
spursConfig.SavePSI = 0;
spursConfig.UseW_Regularization = 0;%%!!!!

%% 1.1 Default
% Angle        = 1,3,5,7,9,11
% Degree       = 3
% Iterations   = 1
% Oversampling = 1

spursConfig.KernelFunctionDegree = 3;
spursConfig.Niterations = 1;
spursConfig.OverGridFactor = 1;

IsEmbedded = 1;
numAngels = [1,3,5,7,9,11];

runScripts(numAngels, ProccesOpt, IsEmbedded, ImageOpt, acquisition_type, spursConfig);

%% 1.2 Spline Degree
% Angle        = 1,3
% Degree       = 5
% Iterations   = 1
% Oversampling = 1

spursConfig.KernelFunctionDegree = 5;
spursConfig.Niterations = 1;
spursConfig.OverGridFactor = 1;

IsEmbedded = 1;
numAngels = [1,3];

runScripts(numAngels, ProccesOpt, IsEmbedded, ImageOpt, acquisition_type, spursConfig);

%% 1.3 Iterations
% Angle        = 1,3
% Degree       = 3,5
% Iterations   = 20
% Oversampling = 1


spursConfig.Niterations = 20;
spursConfig.OverGridFactor = 1;

IsEmbedded = 1;
numAngels = [1,3,5];
for iKernelFunctionDegree = [3,5]
    spursConfig.KernelFunctionDegree = iKernelFunctionDegree;
    runScripts(numAngels, ProccesOpt, IsEmbedded, ImageOpt, acquisition_type, spursConfig);
end


function runScripts(numAngels, ProccesOpt, IsEmbedded, ImageOpt, acquisition_type, spursConfig)

for inumAngels = numAngels
    for j = ProccesOpt
        for ii = IsEmbedded
            for i = ImageOpt
                phantom_type = i;
                data_type = j;
                tic;
                script_reconstruct_images_from_das(acquisition_type , phantom_type , data_type, ii, inumAngels, spursConfig);
                fprintf('data = %d. Embedded = %d. Time = %.2f [sec]\n',data_type, ii, toc/60)
            end
        end
    end
end
end