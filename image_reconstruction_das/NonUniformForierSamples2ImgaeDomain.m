function [imageRecover,imageFFT] = NonUniformForierSamples2ImgaeDomain(scan, fx_mesh, fz_mesh, Gamma, fx, fsx, processType, sumForierDomainFlag, spursConfig)

%% 3. Transform back to image.
%cases:
%       3 - Interpolate the all frequency map. IFFT. Hilbert.
%       4 - Interpolate the signal in baseband. IFFT.
%       5 - Using direct iDFT. Accurate transform. (Time issue)
%       6 - Using NUFFT (Jeff Fessler).
imageFFT = 0;
switch processType
    case 3
        imageRecover = single(InterpFullFreq(Gamma, scan, fx , fz_mesh, fx_mesh));
    case 4
        % Average Gamma Over Same Coordinates
        [fx_mesh, fz_mesh, Gamma] = uniqueCordinates(fx_mesh, fz_mesh, Gamma);
        
        imageRecover = InterpFreqBB(Gamma, scan, fx , fz_mesh, fx_mesh, sumForierDomainFlag);
    case 5
        imageRecover = NUDFT(Gamma, scan, fz_mesh, fx_mesh , fsx);
    case 6
        %         [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz_mesh, fx_mesh , fx, fsx);
        imageRecover = single(InterpNUFFT(Gamma, scan, fz_mesh, fx_mesh , fsx, sumForierDomainFlag));
    case 7
<<<<<<< HEAD
        %         profile on;
        %         [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz_mesh, fx_mesh , fx, fsx);
        [b, Kappa_m,sqrtN] = SPURSInit2(Gamma, scan, fz_mesh, fx_mesh , fx, fsx, sumForierDomainFlag);
=======
%         profile on;
%         [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz_mesh, fx_mesh , fx, fsx);
        [b, Kappa_m,sqrtN, W] = SPURSInit2(Gamma, scan, fz_mesh, fx_mesh , fx, fsx, sumForierDomainFlag);
>>>>>>> 0da018a837fe7c2cedf3360244b33d2dc085cc6a
        if nargin < 9
            BsplineDegree = 3;
            Rho = 1e-3;%%
            Niterations = 1;
            OverGridFactor = 2;
            FilterInImageSpace = 1;
            savePsi = 0;
            UseW_Regularization = 0;
        else
            BsplineDegree = spursConfig.KernelFunctionDegree;
            Rho = spursConfig.Rho;
            Niterations = spursConfig.Niterations;
            OverGridFactor = spursConfig.OverGridFactor;
            FilterInImageSpace = spursConfig.FilterInImageSpace;
            savePsi = spursConfig.SavePSI;
            UseW_Regularization = spursConfig.UseW_Regularization;
        end
        
        SPURS_settings.sqrtN = sqrtN;
        SPURS_settings.KernelFunctionString = 'Bspline';
        SPURS_settings.KernelFunctionDegree = BsplineDegree;
        SPURS_settings.ReusePrecalculatedData = 1;
        SPURS_settings.Rho = Rho;
        SPURS_settings.Niterations = Niterations;
        SPURS_settings.UseW = 0;
        SPURS_settings.ForceGenrateNewPhi = 0;
        SPURS_settings.ForceFactorPsi = 0;
        SPURS_settings.SavePSI = savePsi;
        SPURS_settings.OverGridFactor = OverGridFactor;
        SPURS_settings.alpha = 1;
        SPURS_settings.CalcOptimalAlpha = 1;
        SPURS_settings.FilterInImageSpace = FilterInImageSpace;
        
<<<<<<< HEAD
%         [ OutputImages , b_hat] = SPURS(double(b), double(Kappa_m), SPURS_settings);
        [ OutputImages , b_hat, imageFFT] = SPURS_Tanya(double(b), double(Kappa_m), SPURS_settings);
        
        %         [ OutputImages_real  , b_hat] = SPURS(double(real(b)), double(Kappa_m), SPURS_settings);
        %         [ OutputImages_imag  , b_hat] = SPURS(double(imag(b)), double(Kappa_m), SPURS_settings);
=======
        SPURS_settings.UseW_Regularization = UseW_Regularization;
        SPURS_settings.W_Regularization = W;
        
        [ OutputImages , b_hat] = SPURS(double(b), double(Kappa_m), SPURS_settings);
%         [ OutputImages_real  , b_hat] = SPURS(double(real(b)), double(Kappa_m), SPURS_settings);
%         [ OutputImages_imag  , b_hat] = SPURS(double(imag(b)), double(Kappa_m), SPURS_settings);
>>>>>>> 0da018a837fe7c2cedf3360244b33d2dc085cc6a
        imageRecover = OutputImages(:,:,end);
        %         profile viewer;
    case 8
        %         [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz_mesh, fx_mesh , fx, fsx);
        [imageRecover, imageFFT] = (InterpNUFFT2(Gamma, scan, fz_mesh, fx_mesh , fsx, sumForierDomainFlag));
    case 9
        % Average Gamma Over Same Coordinates
        [fx_mesh, fz_mesh, Gamma] = uniqueCordinates(fx_mesh, fz_mesh, Gamma);
        
        imageRecover = single(InterpLinearSlice(Gamma, scan, fz_mesh, fx_mesh , fsx, sumForierDomainFlag));
end


end
