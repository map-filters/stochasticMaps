% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

function filter_out = post_process(model, filter, options)
% Post_process compute the statistics (mean, median, and 
% standard deviation) of the RMSE, spread and coverage 
% probability over (J-T_BurnIn) assimilation steps.

% declare struct to save results
filter_out = struct;

% extract indices for assimilation cycles
n0 = round(model.t0/(model.dt*model.dt_iter)) + 1;
Acycles = n0:n0+model.J-1;

% extract model information
d  = model.d;
xt = model.xt(:,Acycles);

% set T_BurnIn
if isfield(model,'T_BurnIn')
    T_BurnIn = model.T_BurnIn;
else
    T_BurnIn = 0;
end

% compute root mean square error
xt_mean = filter.mean(:,T_BurnIn+1:end); 
xt_filt = xt(:,T_BurnIn+1:end);
RMSE = sqrt( sum( (xt_mean - xt_filt).^2,1))/sqrt(d);

% compute RMSE statistics
RMSE_median = median(RMSE);
RMSE_mean   = mean(RMSE);
RMSE_std    = std(RMSE);

% compute ensemble spread 
xt_var = filter.var(:,T_BurnIn+1:end);
spread = sqrt(sum(xt_var,1))/sqrt(d);

% compute spread statistics
spread_mean   = mean(spread); 
spread_median = median(spread);
spread_std    = std(spread);

% compute quantile information
qinf = squeeze(filter.quantiles(:,1,T_BurnIn+1:end));
qsup = squeeze(filter.quantiles(:,2,T_BurnIn+1:end));

% compute coverage prob (average over time)
coverage_prob = mean(xt_filt >= qinf & xt_filt <= qsup,1);

% compute coverage prob statistics
covprob_mean   = mean(coverage_prob);
covprob_median = median(coverage_prob);
covprob_std    = std(coverage_prob);

% save results in filter_out
filter_out.RMSE           = RMSE;
filter_out.RMSE_median    = RMSE_median;
filter_out.RMSE_mean      = RMSE_mean;
filter_out.RMSE_std       = RMSE_std;

filter_out.spread         = spread;
filter_out.spread_mean    = spread_mean;
filter_out.spread_median  = spread_median;
filter_out.spread_std     = spread_std;

filter_out.coverage_prob  = coverage_prob;
filter_out.covprob_mean   = covprob_mean;
filter_out.covprob_median = covprob_median;
filter_out.covprob_std    = covprob_std;

% copy options to filter_out
for fn = fieldnames(options)'
   filter_out.(fn{1}) = options.(fn{1});
end

end

% -- END OF FILE --
