 function x = dtft2_adj(X, omega, N1, N2, n_shift, useloop)
%function x = dtft2_adj(X, omega, N1, N2, n_shift, useloop)
% Compute adjoint of 2D DTFT for spectrum X at frequency locations omega
% In
%	X	[M,L]		2D DTFT values
%	omega	[M,2]		frequency locations (radians)
%	n_shift [2,1]		use [0:N-1]-n_shift (default [0 0])
%	useloop			1 to reduce memory use (slower)
% Out
%	x	[N1,N2,L]	signal values
%
% Requires enough memory to store M * (N1*N2) size matrices (for testing)
%
% Copyright 2001-9-17, Jeff Fessler, The University of Michigan

if ~exist('n_shift','var') || isempty(n_shift), n_shift = [0 0]; end
if ~exist('useloop','var') || isempty(useloop), useloop = 0; end

[nn1 nn2] = ndgrid([0:(N1-1)]-n_shift(1), [0:(N2-1)]-n_shift(2));

%
% loop way: slower but less memory
%
if useloop
	M = length(omega);
	x = single(zeros(N1,N2,ncol(X)));		% [N1,N2,M]
	t1 = (1i) * nn1;
	t2 = (1i) * nn2;
	for ii=1:M
		x = x + exp(omega(ii,1)*t1 + omega(ii,2)*t2) * X(ii,:);
        if mod(ii,1e3) == 0
            fprintf('Pass %.2f per\n',ii/M*100)
        end
	end
else
	x = exp(1i*(nn1(:)*omega(:,1)' + nn2(:)*omega(:,2)')) * X; % [N1*N2,L]
	x = reshape(x, [N1 N2 numel(x)/N1/N2]);		% [N1,N2,L]
end
