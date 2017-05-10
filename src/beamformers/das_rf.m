function image = das_rf(scan,dataset,pw_indices)

    %-- Function which implements the conventional Delay And Sum (DAS) beamform technique with apodization in reception
    %-- The corresponding code is dedicated to the reconstrucion of dataset (rawdata) saved in RF format

    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  

    assert(isempty(dataset.modulation_frequency)||dataset.modulation_frequency==0,'The supplied dataset is not RF');

    %-- select the plane waves that will be used in each frame
    if nargin < 3
        pw_indices{1} = 1:dataset.firings;
    end

    %-- define scan based on time axis
    time = (0:(size(dataset.data,1)-1)).'/dataset.sampling_frequency+dataset.initial_time;
    z_axis= time*dataset.c0/2;
    rf_scan = linear_scan(scan.x_axis,z_axis);

    %-- receive apodization
    %-- dynamically expanding receive aperture with hanning apodization
    rx_f_number = 1.75;
    rx_aperture = rf_scan.z/rx_f_number;
    rx_aperture_distance = abs(rf_scan.x*ones(1,dataset.channels)-ones(rf_scan.pixels,1)*dataset.probe_geometry(:,1).');
    receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey25');

    %-- angular apodization -> no apodization
    angular_apodization = ones(rf_scan.pixels,dataset.firings); 

    %-- beamforming loop
    beamformed_data = zeros(rf_scan.pixels,length(pw_indices));
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    wb=waitbar(0,'DAS beamforming');

    for f=1:length(pw_indices)
        waitbar(f/length(pw_indices),wb,sprintf('DAS-RF beamforming %0.0f%%',f/length(pw_indices)*100));
        for pw=pw_indices{f}
            %-- transmit delay
            transmit_delay = rf_scan.z*cos(dataset.angles(pw))+rf_scan.x*sin(dataset.angles(pw));
            for nrx=1:dataset.channels
                %-- receive delay
                receive_delay = sqrt((dataset.probe_geometry(nrx,1)-rf_scan.x).^2+(dataset.probe_geometry(nrx,3)-rf_scan.z).^2);
                %-- total delay
                delay = (transmit_delay+receive_delay)/dataset.c0;
                %-- beamformed data
                beamformed_data(:,f) = beamformed_data(:,f)+angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
            end
            clc;
            disp([num2str(pw),' / ',num2str(length(pw_indices{f}))])           
        end
    end
    close(wb);

    beamformed_data(isnan(beamformed_data))=0;
    
    %% Our Section
    probeIdx = 38; %38
    theta = dataset.angles(probeIdx);
    
    
    fsx = 128/(dataset.probe_geometry(end,1) - dataset.probe_geometry(1,1));
    fx = (-1/2:1/128:1/2-1/128)*fsx;
    f = (-1/2:1/dataset.samples:1/2-1/dataset.samples)*dataset.sampling_frequency;
    [x,y] = meshgrid(fx, f);
    S = abs(fftshift(fft2(dataset.data(:,:,probeIdx))));
    figure
    mesh(x, y, S);
    title('S(f,f_x)','fontsize',24)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    
    
    fIdxDown = floor(numel(f)*(0.5 + 3e6/dataset.sampling_frequency));
    fIdxUp = floor(numel(f)*(0.5 + 7e6/dataset.sampling_frequency));
%     fIdxUp = numel(f) - (fIdxDown - floor(numel(f)/2));
    fIdx = (fIdxDown:fIdxUp);
    fTmp = f(fIdx);
    lambda = dataset.c0./fTmp';%%!!dataset.c0./f
    tmp1 = bsxfun(@minus, fx, sin(theta)./lambda);
    tmp2 = bsxfun(@minus, 1./(lambda.^2), tmp1.^2);
    
    fz = bsxfun(@plus, cos(theta)./lambda, tmp2.^0.5);
    
        tmpNom = bsxfun(@minus, 1./(lambda.^2), fx.^2 ); %%! isin(lambda)
    tmpDenom=  -4*pi*  1./(lambda.^2);
    constS = bsxfun(@times, tmpNom, 1./tmpDenom);
    sRes = S(fIdx,:).*constS;
    figure
    fxRelevant = x(fIdx,:);
    mesh(fxRelevant, fz, sRes);
    title('\Gamma(f,f_x)','fontsize',24)
    xlabel('f_x','fontsize',18)
    ylabel('f_z','fontsize',18)
%     N1 = 2^6;		% # of time samples
%     n = [0:N1-1]'-N1/2;	% time sample locations
%     dtft2_adj(X, omega, N1, N2, n_shift, useloop)
%     xd = dtft2_adj(Stmp, [x(fIdx,:) fz], N1, 1, [N1/2 0]);
    %% Interp1
    delta = 1;
    minInd = floor(min(fz(:)));
    maxInd = ceil(max(fz(:)));
    fzWanted = (minInd: delta :maxInd)';
    
    sResInterp1 = zeros(maxInd, size(fxRelevant,2));
    for i = 1:size(fz,2)
        sResInterp1(minInd:end,i) = interp1(fz(:,i), sRes(:,i), fzWanted);
    end
    fxAxesWanted = fxRelevant(1,:);
    [fxAxesWantedMesh,fzWantedMesh] = meshgrid(fxAxesWanted, (1:maxInd));
    sResInterp1(isnan(sResInterp1)) = 0;
    mesh(fxAxesWantedMesh, fzWantedMesh, sResInterp1);
    shg
    
    % Duplicate
%     sResInterp1Duplicate = [sResInterp1(end:-1:1,:);sResInterp1];
%     [fxAxesWantedMesh,fzWantedMesh] = meshgrid(fxAxesWanted, (-maxInd:maxInd-1));
%     mesh(fxAxesWantedMesh, fzWantedMesh, sResInterp1Duplicate);
    
    imagRecover = ifft2(ifftshift(sResInterp1(4500:8500,:)));
    imagesc(20*log10(abs(imagRecover)/max(abs(imagRecover))))
if(0)
    Nx = floor(scan.Nx);
    Nz = floor(scan.Nz);
    fxImage = 1/scan.dx;
    fzImage = 1/scan.dz;
    [~,maxInd] = max(abs(sRes(:)));
    fzRf = fz(maxInd);
    tic;
    xd = dtft2_adj(sRes(:), [fxRelevant(:)/fxImage*2*pi, (fz(:)-fzRf)/fzImage*2*pi], Nx, Nz,[Nx/2,Nz/2],1);
    timeTake = toc;
    fprintf('DFT %.2f seconds\n',timeTake);
    imagesc(db(abs(xd)))
end
    %%
    %-- reshape
    reshaped_beamformed_data = reshape(beamformed_data,[numel(rf_scan.z_axis) numel(rf_scan.x_axis)  length(pw_indices)]);

    %-- compute envelope
    envelope_beamformed_data = tools.envelope(reshaped_beamformed_data);

    %-- interpolate the requested grid
    resampled_envelope_beamformed_data = zeros(numel(scan.z_axis),numel(scan.x_axis),numel(pw_indices));
    for f=1:length(pw_indices)
        resampled_envelope_beamformed_data(:,:,f) = interp1(rf_scan.z_axis,envelope_beamformed_data(:,:,f),scan.z_axis,'linear',0);
    end

    %-- declare an us_image object to store the beamformed data
    image = us_image('DAS-RF beamforming');
    image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
    image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
    image.algorithm = 'Delay-and-Sum (RF version)';
    image.scan = scan;
    image.number_plane_waves = cellfun('length',pw_indices);
    image.data = resampled_envelope_beamformed_data;
    image.transmit_f_number = 0;
    image.receive_f_number = rx_f_number;
    image.transmit_apodization_window = 'none';
    image.receive_apodization_window = 'Tukey 25%';

end
