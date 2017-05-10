function image = das_iq(scan,dataset,pw_indices)

    %-- Function which implements the conventional Delay And Sum (DAS) beamform technique with apodization in reception
    %-- The corresponding code is dedicated to the reconstrucion of dataset (rawdata) saved in IQ format

    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  

    assert(~isempty(dataset.modulation_frequency)&&dataset.modulation_frequency~=0,'The supplied dataset is not IQ');

    %-- select the plane waves that will be used in each frame
    if nargin < 3
        pw_indices{1} = 1:dataset.firings;
    end

    %-- receive apodization: 
    %-- dynamically expanding receive aperture with Tukey 25% apodization
    rx_f_number = 1.75;
    rx_aperture = scan.z/rx_f_number;
    rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
    receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey25');

    %-- angular apodization -> no apodization
    angular_apodization = ones(scan.pixels,dataset.firings); 

    %-- beamforming loop
    beamformed_data = zeros(scan.pixels,length(pw_indices));
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;
    wb = waitbar(0,'DAS beamforming');

    for f=1:length(pw_indices) 
        waitbar(f/length(pw_indices),wb,sprintf('DAS-IQ beamforming %0.0f%%',f/length(pw_indices)*100));
        for pw=pw_indices{f}
            %-- transmit delay
            transmit_delay = scan.z*cos(dataset.angles(pw))+scan.x*sin(dataset.angles(pw));
            for nrx=1:dataset.channels
                %-- receive delay
                receive_delay = sqrt((dataset.probe_geometry(nrx,1)-scan.x).^2+(dataset.probe_geometry(nrx,3)-scan.z).^2);
                %-- total delay
                delay = (transmit_delay+receive_delay)/dataset.c0;
                %-- phase shift
                phase_shift = exp(1i.*w0*(delay-2*scan.z/dataset.c0));
                %-- beamformed data
                beamformed_data(:,f) = beamformed_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
            end
            clc;
            disp([num2str(pw),' / ',num2str(length(pw_indices{f}))])   
        end
    end
    close(wb);
%     %% Our Section
%     probeIdx = 38; %38
%     theta = dataset.angles(probeIdx);
%     
%     
%     fsx = 128/(dataset.probe_geometry(end,1) - dataset.probe_geometry(1,1));
%     fx = (-1/2:1/128:1/2-1/128)*fsx;
%     f = dataset.modulation_frequency + (-1/2:1/dataset.samples:1/2-1/dataset.samples)*dataset.sampling_frequency;
%     [x,y] = meshgrid(fx, f);
%     S = fftshift(fft2(dataset.data(:,:,probeIdx)));
%     figure
%     mesh(x, y, abs(S));
%     title('S(f,f_x)','fontsize',24)
%     xlabel('f_x','fontsize',18)
%     ylabel('f','fontsize',18)
%     
%     
% %     fIdxDown = floor(numel(f)*(0.5 + 3e6/dataset.sampling_frequency));
% %     fIdxUp = floor(numel(f)*(0.5 + 7e6/dataset.sampling_frequency));
% % %     fIdxUp = numel(f) - (fIdxDown - floor(numel(f)/2));
% %     fIdx = (fIdxDown:fIdxUp);
% %     fTmp = f(fIdx);
%     lambda = dataset.c0./f';
%     tmp1 = bsxfun(@minus, fx, sin(theta)./lambda);
%     tmp2 = bsxfun(@minus, 1./(lambda.^2), tmp1.^2);
%     fz = bsxfun(@plus, cos(theta)./lambda, tmp2.^0.5);
%     
%     tmpNom = bsxfun(@minus, 1./(lambda.^2), fx.^2 ); %%! isin(lambda)
%     tmpDenom=  -4*pi*  1./(lambda.^2);
%     constS = bsxfun(@times, tmpNom, 1./tmpDenom);
%     sRes = S .* constS;
%     
%     figure
%     fxRelevant = x;
%     mesh(fxRelevant, fz, abs(sRes));
%     title('\Gamma(f,f_x)','fontsize',24)
%     xlabel('f_x','fontsize',18)
%     ylabel('f_z','fontsize',18)
% %     N1 = 2^6;		% # of time samples
% %     n = [0:N1-1]'-N1/2;	% time sample locations
% %     dtft2_adj(X, omega, N1, N2, n_shift, useloop)
% %     xd = dtft2_adj(Stmp, [x(fIdx,:) fz], N1, 1, [N1/2 0]);
%     %% Interp1
%     tic;
%     delta = (1/scan.dz ) /scan.Nz;
%     minInd = floor(min(fz(:))/delta);
%     maxInd = ceil(max(fz(:))/delta);
%     minFz = minInd*delta;
%     maxFz = maxInd*delta;
%     fzWanted = (minFz: delta :maxFz)';
%     
%     sResInterp1 = zeros(maxInd + 1, size(fxRelevant,2));
%     for i = 1:size(fz,2)
%         sResInterp1((minInd+1):end,i) = interp1(fz(:,i), sRes(:,i), fzWanted);
%     end
%     sResInterp1(isnan(sResInterp1)) = 0;
%     
%     fxAxesWanted = fxRelevant(1,:);
%     fzAxesWanted = (0:delta:maxFz);
%     [fxAxesWantedMesh,fzWantedMesh] = meshgrid(fxAxesWanted, fzAxesWanted);
%     figure
%     mesh(fxAxesWantedMesh, fzWantedMesh, abs(sResInterp1));
%     title('\Gamma(f,f_x) After Interpolation','fontsize',24)
%     xlabel('f_x','fontsize',18)
%     ylabel('f_z','fontsize',18)
%     shg
%     
%     fzAxesWantedExt = (-maxFz:delta:maxFz);
%     midFzAxes = floor(numel(fzAxesWantedExt)/2);
%     [fxAxesWantedMeshExt,fzWantedMeshExt] = meshgrid(fxAxesWanted, fzAxesWantedExt);
%     sResInterp1Ext = zeros( size(fxAxesWantedMeshExt));
%     % negative reflecation
%     sResInterp1Ext( 1:midFzAxes , 1:end) = conj( sResInterp1(end:-1:2 , end:-1:1) );
%     % positive
%     sResInterp1Ext( (midFzAxes+1):end , 1:end) = sResInterp1(1:end ,1:end);
%     
%     figure
%     mesh(fxAxesWantedMeshExt, fzWantedMeshExt, abs(sResInterp1Ext));
%     title('\Gamma(f,f_x) After Interpolation & Extend','fontsize',24)
%     xlabel('f_x','fontsize',18)
%     ylabel('f_z','fontsize',18)
%     shg
%     
%    
%     imagRecover = ifft2(ifftshift(sResInterp1Ext));
%     
%     %-- compute envelope
%     imagRecoverEnvelope = tools.envelope(abs(imagRecover));
%     figure
%     imagesc(db(imagRecoverEnvelope/max(imagRecoverEnvelope(:))))
%     a = toc; disp(a)
%     
%     
%      %-- interpolate the requested grid
%     resampled_envelope_beamformed_data = zeros(numel(scan.z_axis),numel(scan.x_axis),numel(pw_indices));
%     for f=1:length(pw_indices)
%         resampled_envelope_beamformed_data(:,:,f) = interp1(rf_scan.z_axis,envelope_beamformed_data(:,:,f),scan.z_axis,'linear',0);
%     end
% 
%     %-- declare an us_image object to store the beamformed data
%     image = us_image('DAS-RF beamforming');
%     image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
%     image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
%     image.algorithm = 'Delay-and-Sum (RF version)';
%     image.scan = scan;
%     image.number_plane_waves = cellfun('length',pw_indices);
%     image.data = resampled_envelope_beamformed_data;
%     image.transmit_f_number = 0;
%     image.receive_f_number = rx_f_number;
%     image.transmit_apodization_window = 'none';
%     image.receive_apodization_window = 'Tukey 25%';
% 
%     
% if(0)
%     Nx = floor(scan.Nx);
%     Nz = floor(scan.Nz);
%     tic;
%     xd = dtft2_adj(sRes(:), [2*pi* fxRelevant(:), 2*pi*(fz(:)-mean(fz(:)))], Nx, Nz,[0,0],1);
%     timeTake = toc;
%     fprintf('DFT %.2f seconds\n',timeTake);
%     imagesc(db(abs(xd))),colorbar, shg
% end
% 
% %%
    
    %-- reshape
    envelope_beamformed_data = abs(reshape(beamformed_data,[numel(scan.z_axis) numel(scan.x_axis)  length(pw_indices)]));

    %-- declare an us_image object to store the beamformed data
    image = us_image('DAS-IQ beamforming');
    image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
    image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
    image.algorithm = 'Delay-and-Sum (IQ version)';
    image.scan = scan;
    image.number_plane_waves = cellfun('length',pw_indices);
    image.data = envelope_beamformed_data;
    image.transmit_f_number = 0;
    image.receive_f_number = rx_f_number;
    image.transmit_apodization_window = 'none';
    image.receive_apodization_window = 'Tukey 25%';
    
end
