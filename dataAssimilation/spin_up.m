% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

function model = spin_up(model, M)
% Function runs spin-up phase using perturbed observation EnKF for the
% specified model. The results are saved in the specified data_folder.

% set EnKF options
options = struct;
options.M   = M;
options.rho = 0;

% Set initial condition
model.x0 = (repmat(model.m0,1,M) + sqrtm(model.C0)*randn(model.d, M))';
model.t0 = 0;
model.J  = model.T_SpinUp;

% Run EnKF for T_SpinUp steps
enkf  = enkfPert(model, options);
prior = seq_assimilation(model, @enkf.sample_posterior);

% Check tracking during spin-up
RMSE = sqrt(sum((prior.mean - model.xt(:,1:model.T_SpinUp)).^2,1)/model.d);
fprintf('Spin-up mean RMSE: %f\n', mean(RMSE))

% Set initial condition
model.x0 = prior.xp;
model.t0 = prior.time(end);
model.J  = model.T_Steps;

end

% -- END OF FILE --