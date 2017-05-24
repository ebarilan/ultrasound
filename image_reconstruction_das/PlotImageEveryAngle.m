function PlotImageEveryAngle(imageRecoverI, probeIdxTot,scan)
for idxP = 1:numel(probeIdxTot)
    env = imageRecoverI{idxP};
    im = 20*log10(env./max(env(:)));
    vrange = [-60 0];
    x_lim = [min( scan.x_matrix(:)) max( scan.x_matrix(:))]*1e3;
    z_lim = [min( scan.z_matrix(:)) max( scan.z_matrix(:))]*1e3;
    figure
    
    
    xAx = 1e3*linspace(scan.x_axis(1),scan.x_axis(end),size(im,2));
    zAx = 1e3*linspace(scan.z_axis(1),scan.z_axis(end),size(im,1));
    imagesc(xAx,zAx,im);
    shading flat; colormap gray; caxis(vrange); colorbar; hold on;
    axis equal manual;
    xlabel('x [mm]');
    ylabel('z [mm]');
    set(gca,'YDir','reverse');
    set(gca,'fontsize',16);
    axis([x_lim z_lim]);
    title(sprintf('NUFFT\n Angle %.2f',dataset.angles(probeIdxTot(idxP))),'fontsize',24)
end

end