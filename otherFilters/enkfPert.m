% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef enkfPert
	% Define EnKF object to assimilate observation at time t by generating 
	% analysis samples from the prior samples using the
	% perturbed observation EnKF
	%
	% Notes: - Ensemble KF requires linear observation operator

	properties

		model      % struct containing definition of observations
		options    % struct containing options for EnKF algorithm

	end

	methods
		function EnKF = enkfPert(model, options, varargin)

			% define EnKF object
			p = ImprovedInputParser;
            addRequired(p,'model');
            addRequired(p,'options');
			parse(p, model, options, varargin{:});
			EnKF = passMatchedArgsToProperties(p, EnKF);

		end %endFunction
		% ---------------------------------------------------------------
		% ---------------------------------------------------------------
		function Xpost = sample_posterior(EnKF, Xpr, Yt)

            % extract M and d
            d = EnKF.model.d;
            M = EnKF.options.M;

			% check input dimensions
			if size(Xpr,1) ~= M || size(Xpr,2) ~= d
			    error('Dimension mismatch of samples, M, and d')
			end

			% Load model parameters
			obs_lin  = EnKF.model.obs_lin;
			sigma_y  = EnKF.model.sigma_y;

			% Load EnKF inputs
			if isfield(EnKF.options, 'locRad')
				locCov = gaspari_cohn(d, EnKF.options.locRad, 1);
			else
				locCov = ones(d,d);
			end

			% Determine the number of observations
			n_obs = size(obs_lin,1);

			% Compute sample mean and covariance
			Mpr = mean(Xpr,1);
			Cpr = 1/(M-1)*(Xpr - Mpr)'*(Xpr - Mpr);

			% Inflate and localize covariance
			Cpr = (1 + EnKF.options.rho)*(locCov.*Cpr);

			% Compute Kalman gain
			K = Cpr*obs_lin'/(obs_lin*Cpr*obs_lin' + sigma_y^2*eye(n_obs));

			% Find analysis samples
			Xpost = Xpr' - K*(obs_lin*Xpr' + sigma_y*randn(n_obs,M) - Yt);
			Xpost = Xpost';

		end %endFunction
		% ---------------------------------------------------------------
		% ---------------------------------------------------------------
	end %endMethods
end %endClass