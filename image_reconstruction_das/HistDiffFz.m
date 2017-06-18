function dat = HistDiffFz(fx_mesh, fz_mesh, nAngles)

[C,~,ic] = unique(fx_mesh);
dat = cell(size(C));
for i = 1:numel(C)
    fz_curr = sort(fz_mesh(ic == i));
    diff_fz = fz_curr(2:end) - fz_curr(1:end-1);
    dat_i = [C(i)*ones(size(diff_fz)) , diff_fz];
    dat{i} = dat_i;
end
dat = cat(1,dat{:});
figure
h = histogram2(dat(:,1),dat(:,2),'FaceColor','flat');
h.YBinEdges = [-Inf 0:1:60 Inf];
title(['Histogram of Diff fz. Num of Angles = ', num2str(nAngles)])
xlabel('f_x')
ylabel('\Deltaf_z')
view([-93.1 20.4]);

end