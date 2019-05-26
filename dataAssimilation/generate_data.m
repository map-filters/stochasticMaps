% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

function model = generate_data(model, J)

% Load model inputs
d 		  = model.d;
for_op    = model.for_op;
sigma_x   = model.sigma_x;
dt        = model.dt;
dt_iter   = model.dt_iter;
sampleLik = model.sampleLik;
x0        = model.x0;

% find number of observations in space from linear operator
nobs = length(model.data_idx);

% Declare vectors to store state, observations and time
xt = zeros(d, J);
yt = zeros(nobs, J);
tt = zeros(J, 1);

% Initialize xf and tf
xf = x0;
tf = 0;

% Generate true data
for n=1:J

	% run dynamics and save results
	for i=1:dt_iter
		xf = for_op(xf, tf + dt*(i-1), dt);
	end
	xf = xf + sigma_x*randn(size(xf));
	tf = tf + dt*dt_iter;

	% collect observations
	yt(:,n) = sampleLik(model.obs_lin*xf');

	% Save results in vector
	xt(:,n) = xf;
	tt(n)   = tf;

end

% Save data in model
model.xt     = xt;
model.yt     = yt;
model.tt     = tt;

end

% -- END OF FILE --