% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

function filter_out = run_filter_params(model, method, tune_params, fixed_params, n_proc)
% Run filtering algorithm specified in inferenceAlg for each combination of 
% parameters in tune_params cell and post_process results to identify optimal
% set of parameters. The step is implemented using n_proc parallel processing.

% convert tune_param values and fields to cell
params_cell  = struct2cell(tune_params);
params_field = fieldnames(tune_params);

% number of parameters
n_params = size(params_field,1);

% create grid of parameter combinations
param_grid = combvec(params_cell{:})';
n_tests = size(param_grid,1);

% set parallel for loop
c = parcluster('local');
c.NumWorkers = n_proc;
parpool(n_proc);

% extract parameters
parfor i=1:n_tests
    
    fprintf('Test %d/%d\n', i, n_tests)
    
    % extract parameter values
    param_values = param_grid(i,:);

	% set tuning parameters for algorithm based on combinations in param_grid
    options = struct;
    for j=1:n_params
        options.(params_field{j}) = param_values(j);
    end
    
    % add fixed_params to options
    for fn = fieldnames(fixed_params)'
        options.(fn{1}) = fixed_params.(fn{1});
    end
    
    % run filter and post-process outputs
    algorithm = method(model, options);
    filter = seq_assimilation(model, @algorithm.sample_posterior);
    filter_out(i) = post_process(model, filter, options);
    
end

% delete parallel pool object
delete(gcp('nocreate'))

end

% -- END OF FILE --