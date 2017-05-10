function fx = CreateFreqAxes(N , fsx)
fx = (0:1/N:1-1/N)*2*pi;
fx = wrapToPiB(fftshift(fx))/(2*pi)*fsx;
end