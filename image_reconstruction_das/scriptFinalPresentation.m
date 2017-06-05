ImageOpt = 1;%[1,2];
ProccesOpt = 7;%[4, 6]; %[1, 4, 6];
numAngels = 1;
IsEmbedded = 1;
acquisition_type = 1;

spursConfig.KernelFunctionDegree = 3;
spursConfig.Rho = 1e-3;
spursConfig.Niterations = 1;
spursConfig.OverGridFactor = 1;%2;
spursConfig.FilterInImageSpace = 1;
spursConfig.SavePSI = 0;


for i = ImageOpt
    for ii = IsEmbedded
        for inumAngels = numAngels
            for j = ProccesOpt
                phantom_type = i;
                data_type = j;
                tic;
                script_reconstruct_images_from_das(acquisition_type , phantom_type , data_type, ii, inumAngels, spursConfig);
                fprintf('data = %d. Embedded = %d. Time = %.2f [sec]\n',data_type, ii, toc/60)
            end
        end
    end
end