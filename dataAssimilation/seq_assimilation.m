% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

function pf = seq_assimilation(model, analysis_alg)
% SEQ_ASSIMILATION: Function assimilates observations sequentially
% using the algorithm specified with the function handle in analysis_alg
%
% Example analysis_alg: enkf_pert, enkf_pert_seqobs, StochasticMapFilter

% Load model inputs
for_op  = model.for_op;
sigma_x = model.sigma_x;
d       = model.d;
J       = model.J;
dt      = model.dt;
x0      = model.x0;
t0      = model.t0;
dt_iter = model.dt_iter;

% Declare vectors to store results
mean_t 		= zeros(d,J);
var_t  		= zeros(d,J);
cov_t  		= zeros(d,d,J);
quantiles 	= zeros(d,2,J); %store [0.025, .975] quantiles
time_t 		= zeros(J,1);

% Set xf and tf based on initial ensemble of particles
xf = x0;
tf = t0;

% set indices for assimilation cycles
n0 = round(t0/(model.dt*model.dt_iter)) + 1;
Acycles = n0:n0+J-1;

% Run particle filter
for n=1:length(Acycles)

	% -------------------------------------------- %
	% FORECAST 									   %
	% -------------------------------------------- %
	
	% propagate samples over dt_iter time_steps
	for i=1:dt_iter
		xf = for_op(xf, tf + dt*(i-1), dt);
	end
	xf = xf + sigma_x*randn(size(xf));
	tf = tf + dt*dt_iter;
	
	% -------------------------------------------- %
	% ASSIMILATION 								   %
	% -------------------------------------------- %
    
    % extract data at assimilation step n
    Yt = model.yt(:,Acycles(n));

	% generate posterior samples
	xp = analysis_alg(xf, Yt);
        
	% -------------------------------------------- %
	% SAVE RESULTS 								   %
	% -------------------------------------------- %

	% save time step
	time_t(n) = tf;

	% compute and save statistics
	mean_t(:,n)      = mean(xp,1);
	var_t(:,n)       = var(xp,[],1);
	cov_t(:,:,n)     = cov(xp);
	quantiles(:,:,n) = quantile(xp,[.025,.975],1)';

	% set next forecast samples
	xf = xp;

end

% Return results
pf = struct;
pf.xp        = xp;
pf.time      = time_t;
pf.mean      = mean_t;
pf.var       = var_t;
pf.cov       = cov_t;
pf.quantiles = quantiles;

% -- END OF FILE --