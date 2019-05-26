% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

close all; clc
sd = 10; rng(sd);

%% ----------------------------
%% SETUP DA PROBLEM
%% ----------------------------

% Define model parameters
T_Steps   = 2000;
T_BurnIn  = 1000;
T_SpinUp  = 2000;

d         = 3;
beta      = 8/3;
rho       = 28;
sigma     = 10;

sigma_x   = 1e-2;
sigma_y   = 2;

dt        = 0.05;
dt_iter   = 2;
ds_obs    = 1;

loglik    = @(x,y) log_scalar_Gaussian(x,y,sigma_y);
sampleLik = @(x) x(:) + sigma_y*randn(length(x),1);

% set initial condition for data generation & spin-up
m0 = zeros(d,1);
C0 = eye(d);

% Setup forward operator
for_op = @(xt, t, dt) rk4(@(t,y) lorenz63(t,y,sigma,rho,beta), xt, t, dt);

% Setup observation operator
obs_s  = 1:ds_obs:d;
nobs_s = length(obs_s);
H = eye(d); H = H(obs_s,:);
obs_op = @(x) H*x;

% Save model parameters
model = struct;
model.d         = d;
model.dt        = dt;
model.dt_iter   = dt_iter;
model.sigma_x   = sigma_x;
model.sigma_y   = sigma_y;
model.loglik    = loglik;
model.sampleLik = sampleLik;
model.data_idx  = obs_s;
model.m0        = m0;
model.C0        = C0;

% save number of steps
model.T_BurnIn  = T_BurnIn;
model.T_Steps   = T_Steps;
model.T_SpinUp  = T_SpinUp;

% Save operators
model.for_op    = for_op;
model.ds_obs    = ds_obs;
model.obs_op    = obs_op;
model.obs_lin   = H;

% Save simulation parameters
model.seed      = sd;

fprintf('Setup DA problem\n')

%% ----------------------------
%% GENERATE DATA
%% ----------------------------

% set initial condition
model.x0 = (m0 + sqrtm(C0)*randn(d,1))';

% run dynamics and generate data
model = generate_data(model, T_SpinUp + T_Steps);

fprintf('Generated data\n')

% -- END OF FILE --