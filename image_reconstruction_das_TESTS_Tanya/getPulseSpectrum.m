function [pulseSpectrum] = getPulseSpectrum(f,f0,fBW)

B = f0*fBW; 
a = (pi*B/2)^2/log(2); 
% transducerFreqResponseTheor =
% 0.5*sqrt(pi/a)*(exp(-pi^2/a*(f-f0).^2)+exp(-pi^2/a*(f+f0).^2)); % two-sided
transducerFreqResponseTheor = 0.5*sqrt(pi/a)*(exp(-pi^2/a*(f-f0).^2)); % one-sided
pulseSpectrum = transducerFreqResponseTheor/max(transducerFreqResponseTheor);
figure; plot(f/1e6,abs(pulseSpectrum));
string = ['Gaussian pulse, f_0 = ', num2str(f0/1e6), ' MHz ']; title(string);
xlabel('f [MHz]'); 


