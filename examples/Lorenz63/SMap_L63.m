% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

sd = 10; rng(sd);

%% ----------------------------
%% RUN STOCHASTIC MAP FILTER
%% ----------------------------

for orderTM = orderTM_vect
    for M = M_vect

        fprintf('order %d, M = %d\n', orderTM, M);

        % load spin-up data
        load(['data/spinup_M' num2str(M)],'model')

        % set fixed options
        options = struct;
        options.M                = M;
        options.distMat          = metric_lorenz(model.d);
        options.order_all        = orderTM;
        options.nonId_radius     = model.d;
        options.offdiag_rad      = model.d;
        
        % set tuning options
        tune_options = struct;
        tune_options.rho         = rho_all;

        % run filter for all parameters
        filter = run_filter_params(model, @StochasticMapFilter, ...
            tune_options, options, n_proc);

        % save workspace
        save(['data/SM_order' num2str(orderTM) '_M' num2str(M)], 'filter');

    end
end

% -- END OF FILE --