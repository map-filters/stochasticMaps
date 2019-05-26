% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

clear; close all; clc
sd = 10; rng(sd);

% add functions to current path
addpath('../../dataAssimilation')
addpath('../../methods')
addpath('../../otherFilters')
addpath('../../tools')

% define number of samples
M_vect = [10,20,40,60,80,100,200,400,600];

% define tuning parameter for inflation factor
rho_all     = 0:0.05:0.5;
locRad_all  = 0:20;
nonId_all   = 0:10;
offdiag_all = 0:10;

% define orders for stochastic maps
orderTM_vect = [1,2,3];

% define number of processors
n_proc = 10;

% setup L96 problem
L96_problem_setup;
model_init = model;

% Run spin-up and save data
for M = M_vect
	model = spin_up(model_init, M);
    save(['data/spinup_M' num2str(M)],'model')
end
disp('Spin-up complete!')

% Run EnKF
Seq_EnKF
disp('EnKF complete!')

% Run SMaps
SMap_L96
disp('Stochastic Map complete!')

% run post-processing
generate_plots

% -- END OF FILE --
