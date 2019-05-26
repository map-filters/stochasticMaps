% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef enkfSampleObs
	% Define EnKF object to Assimilate observation at time t by 
	% generating analysis samples from the prior samples using 
	% the sequential perturbed observation EnKF. In this EnKF,
	% the inflation factor is applied to the prior samples before
	% sampling from the likelihood and the covariances are computed
	% using inflated prior and (not-inflated) likelihood samples.
	%
	% Notes: - Ensemble KF requires linear observation operator

	properties

		model      % struct containing definition of observations
		options    % struct containing options for EnKF algorithm

	end

	methods
		function EnKF = enkfSampleObs(model, options, varargin)

			% define EnKF object
			p = ImprovedInputParser;
            addRequired(p,'model');
            addRequired(p,'options');
			parse(p, model, options, varargin{:});
			EnKF = passMatchedArgsToProperties(p, EnKF);

		end %endFunction
		% -----------------------------------------------------------------
		% -----------------------------------------------------------------
		function Xpost = sample_posterior(EnKF, Xpr, Yt)

			% Load model inputs
			obs_lin   = EnKF.model.obs_lin;
			sampleLik = EnKF.model.sampleLik;

			% Load EnKF inputs
			rho       = EnKF.options.rho;

			% extract subset of data and data_idx
			data_idx  = EnKF.model.data_idx;

			% Determine the number of observations
			n_obs = size(obs_lin,1);

			% set gaspari_cohn matrix
			if isfield(EnKF.options, 'locRad')
				GC_mat = gaspari_cohn(EnKF.model.d, EnKF.options.locRad, 1);
			else
				GC_mat = ones(EnKF.model.d,EnKF.model.d);
			end

			% run assimilation sequentially
			xf = Xpr;
			for i=1:n_obs

				% extract local data and local observation operator
				obs_lin_local = obs_lin(i,:);
				data_local = Yt(i,:);

				% sample from conditional likelihood using inflated particles
				mean_x = mean(xf,1);
				x_infl = sqrt(1 + rho)*(xf - mean_x) + mean_x;
				y_sample = sampleLik(obs_lin_local*x_infl');
				
				% join x,y samples
				xy_samples = [y_sample, x_infl];
				mean_xy    = mean(xy_samples,1);
				delta_xy   = xy_samples - mean_xy;
				
				% estimate covariances
				CovAll = 1/(EnKF.options.M-1)*(delta_xy'*delta_xy);
				CovXY  = CovAll(2:end,1);
				VarY   = CovAll(1,1);

				% localize covariance and compute Kalman gain
				CovXY = CovXY.*GC_mat(:,data_idx(i));
				K = CovXY/VarY;
			    
				% sample from conditional likelihood using un-inflated particles
				y_sample = sampleLik(obs_lin_local*xf');

				% Find analysis samples
				xf = (xf' - K*(y_sample' - data_local))'; 

			end

			% set Xpost
			Xpost = xf;

		end %endFunction
		% -----------------------------------------------------------------
		% -----------------------------------------------------------------
	end %endMethods
end %endClass