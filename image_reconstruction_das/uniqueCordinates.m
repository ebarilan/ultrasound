function [fx_mesh_unique, fz_mesh_unique, Gamma_unique] = uniqueCordinates(fx_mesh, fz_mesh, Gamma)
% Average Gamma Over Same Coordinates
[ux,~,idx] = unique([fx_mesh(:), fz_mesh(:)],'rows');
GammaMean = accumarray(idx,Gamma(:),[],@mean);
fx_mesh_unique = ux(:,1);  fz_mesh_unique = ux(:,2);  Gamma_unique = GammaMean;