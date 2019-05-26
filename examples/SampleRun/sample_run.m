% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

% Sample run for Lorenz 63 model over 2000 assimilation steps
% using the Stochastic Map algorithm

clear; close all; clc
sd = 10; rng(sd);

% add paths
addpath('../../dataAssimilation')
addpath('../../methods')
addpath('../../otherFilters')
addpath('../../tools')

% define parameters 
orderTM = 2;    % map order
M       = 200;  % number of samples
rho     = 0.1;  % inflation factor

%% ----------------------------
%% RUN TMAP FILTER
%% ----------------------------

% run spin-up for Lorenz-63
L63_problem_setup
model = spin_up(model, M);

% set fixed options
options = struct;
options.M            = M;
options.distMat      = metric_lorenz(model.d);
options.order_all    = orderTM;
options.nonId_radius = model.d;
options.offdiag_rad  = model.d;
options.rho          = rho;

% run filter
SM = StochasticMapFilter(model, options);
filter = seq_assimilation(model, @SM.sample_posterior);

% plot filter mean over time
figure('position',[0,0,1200,1200])
for i = 1:model.d
    subplot(model.d,1,i)
    hold on
    plot(1:model.J, model.xt(i,model.T_SpinUp+1:end), '-k')
    plot(1:model.J, filter.mean(i,:), '--b')
    xlabel('Time $t$')
    ylabel(['$x_{' num2str(i) '}$'])
    xlim([model.T_BurnIn, model.J])
    hold off
end

% -- END OF FILE --