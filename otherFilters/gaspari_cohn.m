% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

function G = gaspari_cohn(d, L, periodic)
% GASPARI_COHN: Function returns the Gaspari-Cohn tapering
% matrix used for localizing covariance matrices in EnKF based
% on a localization radius L

G = zeros(d,d);
for i=1:d

	% check for periodic domain
	if periodic == 1
		r_dom  = min([abs((1:d) - i); abs((-d:-1) - i); abs((d+1:2*d) - i)]);
	else
		r_dom  = abs((1:d) - i);
	end

	% evaluate GC_row
	G(i,:) = GC_row(r_dom/L);

end

end

function GC_eval = GC_row(ri)
% GC_ROW: Evaluate the decay based on the Gaspari-Cohn function

	% declare vector to store evaluations of GC function
	GC_eval = zeros(length(ri),1);

	% Find indices for ranges of r
	idx_l1  = (0 <= ri) & (ri < 1);
	idx_1l2 = (1 <= ri) & (ri < 2);
	idx_g2  = (ri >= 2);

	% Find entries of ri in each range
	r_l1  = ri(idx_l1);
	r_1l2 = ri(idx_1l2);
	r_g2  = ri(idx_g2);

	% Compute entries of G
	G_l1  = 1 - 5/3*r_l1.^2 + 5/8*r_l1.^3 + 0.5*r_l1.^4 - 1/4*r_l1.^5;
	G_1l2 = 4 - 5*r_1l2 + 5/3*r_1l2.^2 + 5/8*r_1l2.^3 - 0.5*r_1l2.^4 + 1/12*r_1l2.^5 - 2/3./r_1l2;
	G_g2  = zeros(length(r_g2),1);

	% Assemble GC_eval
	GC_eval(idx_l1)  = G_l1;
	GC_eval(idx_1l2) = G_1l2;
	GC_eval(idx_g2)  = G_g2;

end

% -- END OF FILE --