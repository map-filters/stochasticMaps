% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

sd = 10; rng(sd);

%% ----------------------------
%% RUN SEQUENTIAL EnKF FILTER
%% ----------------------------

for M = M_vect

    fprintf('M = %d\n', M);

    % load data
    load(['data/spinup_M' num2str(M)],'model')

    % define fixed and tuning options
    options = struct;
    options.M = M;

    tune_options = struct;
    tune_options.locRad = locRad_all;
    tune_options.rho    = rho_all;

    % run filter for all parameters
    filter = run_filter_params(model, @enkfSampleObs, tune_options, options, n_proc);

    % Save workspace
    save(['data/enkf_pert_M' num2str(M)], 'filter');

end

% -- END OF FILE --
