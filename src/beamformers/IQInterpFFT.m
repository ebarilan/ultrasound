function [image,imageFFT] = IQInterpFFT(scan,dataset,pw_indices, processType, isEmbedded, spursConfig)
%% 1. Transform the signal to Fourier domain [f,fx]
probeIdxTot = pw_indices{1};
nAngles = numel(probeIdxTot);

imageFFT = 0;
if isEmbedded
    fx_mesh_ceil = cell(nAngles, 1); fz_mesh_ceil = cell(nAngles, 1); Gamma_ceil = cell(nAngles, 1);
    for i = 1:nAngles
        probeIdx = probeIdxTot(i);
        
        [fx_mesh_tmp, fz_mesh_tmp, Gamma_tmp, fx, fsx] = Sampels2FourierDomain(dataset,probeIdx);
%         [fx_mesh_tmp, fz_mesh_tmp, Gamma_tmp, fx, fsx] = Sampels2FourierDomainTanya(dataset,probeIdx);

        fx_mesh_ceil{i} = fx_mesh_tmp(:); fz_mesh_ceil{i} = fz_mesh_tmp(:); Gamma_ceil{i} = Gamma_tmp(:);
    end
    
    fx_mesh = cat(1,fx_mesh_ceil{:}); fz_mesh = cat(1,fz_mesh_ceil{:}); Gamma =cat(1,Gamma_ceil{:});
    
    
    marginOnly = 0;%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!
    
    if marginOnly
        % Embbeding only the frequency margin
        midIndAngle = ceil(nAngles/2);
        minFx = min(fx_mesh_ceil{midIndAngle});
        maxFx = max(fx_mesh_ceil{midIndAngle});
        cond = @(x) ((x < minFx) | (x > maxFx));
        fx_mesh_margin_only = cell(nAngles,1); fz_mesh_margin_only = cell(nAngles,1); Gamma_margin_only = cell(nAngles,1);
        for i = 1:nAngles
            if i == midIndAngle
                fx_mesh_margin_only{i} = fx_mesh_ceil{i}; fz_mesh_margin_only{i} = fz_mesh_ceil{i}; Gamma_margin_only{i} = Gamma_ceil{i};
%                 fx_mesh_margin_only{i} = repmat(fx_mesh_margin_only{i},floor(nAngles/2/10),1); fz_mesh_margin_only{i} = repmat(fz_mesh_margin_only{i},floor(nAngles/2/10),1); Gamma_margin_only{i} = repmat(Gamma_margin_only{i},floor(nAngles/2/10),1);
                continue; 
            end
            fx_mesh_margin_only{i} = fx_mesh_ceil{i}(cond(fx_mesh_ceil{i})); fz_mesh_margin_only{i} = fz_mesh_ceil{i}(cond(fx_mesh_ceil{i})); Gamma_margin_only{i} = Gamma_ceil{i}(cond(fx_mesh_ceil{i}));
        end
        
        fx_mesh = cat(1,fx_mesh_margin_only{:}); fz_mesh = cat(1,fz_mesh_margin_only{:}); Gamma =cat(1,Gamma_margin_only{:});        
    end
    
    HistDiffFz(fx_mesh, fz_mesh, nAngles);
    
    [imageRecover,imageFFT] = NonUniformForierSamples2ImgaeDomain(scan, fx_mesh, fz_mesh, Gamma, fx, fsx, processType, isEmbedded, spursConfig); 
    
else
    imageRecoverI = cell(nAngles, 1);
    for i = 1:nAngles
        probeIdx = probeIdxTot(i);
        [fx_mesh, fz_mesh, Gamma, fx, fsx] = Sampels2FourierDomain(dataset,probeIdx);
%         [fx_mesh, fz_mesh, Gamma, fx, fsx] = Sampels2FourierDomainTanya(dataset,probeIdx);
        imageRecoverI{i} = NonUniformForierSamples2ImgaeDomain(scan, fx_mesh, fz_mesh, Gamma, fx, fsx, processType, isEmbedded, spursConfig); 
    end
    if(0) % Print Plot Result Of Every Angle 
        PlotImageEveryAngle(imageRecoverI, probeIdxTot,scan,dataset); 
    end
    
    imageRecoverI = cat(3,imageRecoverI{:});
    imageRecover = mean(imageRecoverI, 3);
end

%% 4. Image Interpolation
if(0)
    zAxesFinal = single((0:size(imageRecover,1)-1)/(size(imageRecover,1)-1)*scan.z_axis(end));
    % dxCalc = 1 /(2*(fx(end) - fx(1)));
    % if processType==7
    xAxesFinal = ((0:size(imageRecover,2)-1) - size(imageRecover,2)/2)/size(imageRecover,2)*2*scan.x_axis(end);
    % else
    %     xAxesFinal = ((0:size(imageRecover,2)-1) - size(imageRecover,2)/2) *dxCalc;
    % end
    
    % figure
    % imagesc(xAxesFinal , zAxesFinal , db(imageRecover/max(imageRecover(:))))
    % title('Final Image Interp1','fontsize',24)
    % xlabel('X','fontsize',18)
    % ylabel('Z','fontsize',18)
    
    %-- interpolate the requested grid
    [xAxesFinalMesh,zAxesFinalMesh] = meshgrid(single(xAxesFinal),single(zAxesFinal));
    [x_axisMesh , z_axisMesh] = meshgrid(scan.x_axis , scan.z_axis);
    
    resampled_envelope_beamformed_data = zeros(numel(scan.z_axis),numel(scan.x_axis),numel(pw_indices));
    for f=1:length(pw_indices)
        resampled_envelope_beamformed_data(:,:,f) = interp2(xAxesFinalMesh,zAxesFinalMesh, imageRecover(:,:,f),x_axisMesh , z_axisMesh);
    end
    
    % No_resampled_envelope_beamformed_data = imageRecover(((size(imageRecover,1)-scan.Nz+1):end),:);
else
    resampled_envelope_beamformed_data = imageRecover; %%%%!!!!!
end
%-- declare an us_image object to store the beamformed data
image = us_image('DAS-RF beamforming');
image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
image.algorithm = 'Delay-and-Sum (RF version)';
image.scan = scan;
image.number_plane_waves = cellfun('length',pw_indices);
image.data = resampled_envelope_beamformed_data; % image.data = resampled_envelope_beamformed_data;
image.transmit_f_number = 0;
% image.receive_f_number = 0; %%rx_f_number;
image.transmit_apodization_window = 'none';
image.receive_apodization_window = 'Tukey 25%';


titleName = [];

if isEmbedded
    titleName = 'Embedded. ';
end
   
switch processType
    case 4
        titleName = [titleName,'Linear Interpolation'];
    case 6
        titleName = [titleName,'NUFFT'];
    case 7
        titleName = [titleName,'SPURS'];
    case 8
        titleName = [titleName,'NUFFT 1D'];
    case 9
        titleName = [titleName,'Linear Interpolation 1D'];
end
image.name = titleName;

end