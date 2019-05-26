% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

%% ----------------------------
%% SETUP DA PROBLEM
%% ----------------------------

% Define model parameters
T_Steps   = 4000;
T_BurnIn  = 2000;
T_SpinUp  = 2000;

d  	  = 40;
sigma_x   = 0;
sigma_y   = sqrt(0.5);

dt        = 0.01;
dt_iter   = 40;
ds_obs    = 2;

% set initial condition for data generation & spin-up
m0 = zeros(d,1);
C0 = eye(d);

% Setup forward operator
F = 8;
for_op = @(xt, t, dt) rk4(@(t,y) lorenz96(t,y,F), xt, t, dt);

% Setup observation operator
obs_s = 1:ds_obs:d;
n_obs = length(obs_s);
H = eye(d); H = H(obs_s,:);
obs_op = @(x) H*x;

% define likelihood
loglik    = @(x,y) log_scalar_Gaussian(x,y,sigma_y);
sampleLik = @(x) x(:) + sigma_y*randn(length(x),1);

% Save model parameters
model = struct;
model.d          = d;
model.dt         = dt;
model.sigma_x    = sigma_x;
model.sigma_y    = sigma_y;
model.loglik     = loglik;
model.sampleLik  = sampleLik;
model.data_idx   = obs_s;
model.m0         = m0;
model.C0         = C0;

% save number of steps
model.T_BurnIn   = T_BurnIn;
model.T_Steps    = T_Steps;
model.T_SpinUp   = T_SpinUp;

% Save operators
model.for_op     = for_op;
model.dt_iter    = dt_iter;
model.ds_obs     = ds_obs;
model.obs_op     = obs_op;
model.obs_lin    = H;

% Save simulation parameters
model.seed       = sd;

fprintf('Setup DA problem\n')

%% ----------------------------
%% GENERATE DATA
%% ----------------------------

% set initial condition
model.x0 = (m0 + sqrtm(C0)*randn(d,1))';

% run dynamics and generate data
model = generate_data(model, T_SpinUp + T_Steps);

fprintf('Generated Data\n')

% -- END OF FILE --
