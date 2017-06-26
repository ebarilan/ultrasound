function [fx_mesh, fz_mesh, Gamma, fx, fsx] = Sampels2FourierDomain(dataset,probeIdx)
%% 1. Transform the signal to Fourier domain [f,fx]

theta = single(dataset.angles(probeIdx));

Nchannels = dataset.channels;

fsx = Nchannels/single((dataset.probe_geometry(end,1) - dataset.probe_geometry(1,1))); % fx = 1/dx

fx =  CreateFreqAxes(Nchannels , fsx); % Create freq axes [-pi,pi]
f = dataset.modulation_frequency + CreateFreqAxes(dataset.samples , dataset.sampling_frequency); % f = modolation + IQ
f = single(f);
[fx_mesh,f_mesh] = meshgrid(fx, f);

S_no_wrap = fftshift(fft2(dataset.data(:,:,probeIdx)));
fx_mesh_no_wrap = fx_mesh;

if(0)
    figure
    scatter(fx_mesh(:), f_mesh(:), [], abs(S_no_wrap(:)), 'filled');
    title({'S(f,f_x)'; sprintf('Angle %.2f',theta)},'fontsize',16)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    xlim([-4e3,4e3])
end

lambda = single(dataset.c0)./f';
S_f_x = fftshift( fft(dataset.data(:,:,probeIdx),[],1) ,1);

% Hamming Window
windowOption = 1;
if windowOption
    S_F = fftshift( fft(S_f_x,[],2) ,2);
    h = hamming(numel(fx)).';
    S_F_Windowed = bsxfun(@times, S_F, h);
    S_f_x_Windowed = ifft( ifftshift( S_F_Windowed,2),[],2);
    S_f_x = S_f_x_Windowed;
end
H_raw = repmat(h,size(S_F,1), 1);
H_f_x = ifft( ifftshift( H_raw,2),[],2);

fx0 = sin(theta)./lambda;
xGeo = linspace( dataset.probe_geometry(1,1), dataset.probe_geometry(end,1), size(dataset.probe_geometry,1));
expPhase = exp(1i*2*pi * fx0 * xGeo );
S = fftshift( fft(S_f_x.*expPhase,[],2) , 2);
S = single(S);

H_S = fftshift( fft(H_f_x.*expPhase,[],2) , 2);

if theta <= 0
    unwrapLine = fsx/2 + fx0;
    indWrap = (fx_mesh > unwrapLine);
    fx_mesh(indWrap) = fx_mesh(indWrap) - fsx;
else
    unwrapLine = -fsx/2 + fx0;
    indWrap = (fx_mesh < unwrapLine);
    fx_mesh(indWrap) = fx_mesh(indWrap) + fsx;
end

if(0)
    figure
%     subplot(1,2,1)
    scatter(fx_mesh_no_wrap(:), f_mesh(:), [], abs(S(:)), 'filled');
    title({'S(f,f_x-f_{x0}) Shifted'; sprintf('Angle %.2f',theta)},'fontsize',16)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    xlim([-4e3,4e3])
    hold on;
    
    scatter(unwrapLine, f,[], 'd', 'filled')
    hold off;
    
%     subplot(1,2,2)
    figure
    scatter(fx_mesh(:), f_mesh(:),[], abs(S(:)), 'filled');
    title({'S(f,f_x-f_{x0}) Shifted. Spectral Unwrapping'; sprintf('Angle %.2f',theta)},'fontsize',16)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    xlim([-4e3,4e3])
    hold on;
    
    scatter(unwrapLine, f,[], 'd', 'filled')
    hold off;
end
if(0)
%% Hamming Window
windowOption = 0;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% R = rotx(double(theta)*360/(2*pi))
h2d = repmat(hamming(size(S,2)).',size(S,1),1);
R * h2
h2dRotate = fftshift( fft( ifft(ifftshift(h2d,2),[],2).*expPhase,[],2) , 2);
if windowOption
    h = hamming(size(S,2)).';
    [~,indSort] = sort(fx_mesh,2);
    S_Window = zeros(size(S));
    [~,indSortReverse] = sort(indSort,2);
    for i = 1:size(S,1)
        S_Window(i,:) = S(i,:) .* h(indSortReverse(i,:));
    end
 
end
end
%% 2. Transform from [f,fx] -> [fz,fx]

%2.1 fz Calculation

tmp1 = fx_mesh - f_mesh * sin(theta) / dataset.c0;
tmp2 = ((f_mesh / dataset.c0).^2 - tmp1.^2).^0.5;
fz_mesh = f_mesh * cos(theta) / dataset.c0 + tmp2;

fxminusfx0 = bsxfun(@minus, fx_mesh ,fx0);
tmpNom = ( (f_mesh / dataset.c0).^2  -  fxminusfx0.^2 ).^0.5;
tmpDenom = -4*pi* (f_mesh / dataset.c0).^2;
constS = tmpNom ./ tmpDenom;
Gamma = S .* constS;

H_Gamma = H_S .* constS;

if(0)
    figure
    scatter(fx_mesh(:), fz_mesh(:),[], abs(Gamma(:)),'filled')
    title({'\Gamma(f_x,f_z)'; sprintf('Angle %.2f',theta)},'fontsize',16,'Interpreter','tex')
    xlabel('f_x','fontsize',18)
    ylabel('f_z','fontsize',18)
    grid on
    xlim([-4e3,4e3])
    
end

end