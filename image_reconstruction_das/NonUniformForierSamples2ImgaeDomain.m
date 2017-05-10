function imageRecover = NonUniformForierSamples2ImgaeDomain(scan, fx_mesh, fz_mesh, Gamma, fx, fsx, processType)

%% 3. Transform back to image.
%cases:
%       3 - Interpolate the all frequency map. IFFT. Hilbert.
%       4 - Interpolate the signal in baseband. IFFT.
%       5 - Using direct iDFT. Accurate transform. (Time issue)
%       6 - Using NUFFT (Jeff Fessler).

switch processType
    case 3
        imageRecover = single(InterpFullFreq(Gamma, scan, fx , fz_mesh, fx_mesh));
    case 4
        imageRecover = InterpFreqBB(Gamma, scan, fx , fz_mesh, fx_mesh);
    case 5
        imageRecover = NUDFT(Gamma, scan, fz_mesh, fx_mesh , fsx);
    case 6
        [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz_mesh, fx_mesh , fx, fsx);
        imageRecover = single(InterpNUFFT(Gamma, scan, fz_mesh, fx_mesh , fsx,Kappa_m,sqrtN));
    case 7
%         profile on;
        [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz_mesh, fx_mesh , fx, fsx);
        BsplineDegree = 3;
        Rho = 1e-3;%%
        Niterations = 2;
        OverGridFactor = 2;
        FilterInImageSpace = 1;
        
        SPURS_settings.sqrtN = sqrtN;
        SPURS_settings.KernelFunctionString = 'Bspline';
        SPURS_settings.KernelFunctionDegree = BsplineDegree;
        SPURS_settings.ReusePrecalculatedData = 1;
        SPURS_settings.Rho = Rho;
        SPURS_settings.Niterations = Niterations;
        SPURS_settings.UseW = 0;
        SPURS_settings.ForceGenrateNewPhi = 0;
        SPURS_settings.ForceFactorPsi = 0;
        SPURS_settings.SavePSI = 0;
        SPURS_settings.OverGridFactor = OverGridFactor;
        SPURS_settings.alpha = 1;
        SPURS_settings.CalcOptimalAlpha = 1;
        SPURS_settings.FilterInImageSpace = FilterInImageSpace;
        
        [ OutputImages , b_hat] = SPURS(double(b), double(Kappa_m), SPURS_settings);
%         [ OutputImages_real  , b_hat] = SPURS(double(real(b)), double(Kappa_m), SPURS_settings);
%         [ OutputImages_imag  , b_hat] = SPURS(double(imag(b)), double(Kappa_m), SPURS_settings);
        imageRecover = OutputImages(:,:,end);
%         profile viewer;
end


end
