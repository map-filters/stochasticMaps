# Stochastic Maps for Nonlinear Filtering

## What is the Stochastic Map (SM) Algorithm?

The SM Algorithm performs sequential Bayesian inference in a non-Gaussian state-space models with intractable transition kernels. At each assimilation step, the algorithm estimates a potentially nonlinear transformation and uses it to push forward prior (forecast) samples to posterior (analysis) samples. This transformation is based on a coupling between the prior and the posterior distributions at each assimilation step. More information on the algorithm can be found in the [preprint](https://arxiv.org/pdf/1907.00389.pdf).

## Authors

Alessio Spantini (MIT), Ricardo Baptista (MIT), and Youssef Marzouk (MIT)

E-mails: <spantini@mit.edu>, <rsb@mit.edu>, and <ymarz@mit.edu>

## Installation

The SM algorithm is implemented using MATLAB classes and does not require additional packages. The classes for Bayesian Inference cam be found in the `methods` folder. Scripts for data assimilation functions and other filtering algorithms can be found in the `dataAssimilation` and `otherFilters` folders, respectively.

## Example running SM on the Lorenz 63 problem

We provide an example of running the SM algorithm on the chaotic Lorenz 63 problem (defined in [Lorenz 63](https://journals.ametsoc.org/doi/pdf/10.1175/1520-0469%281963%29020%3C0130%3ADNF%3E2.0.CO%3B2)). We consider a configuration with inter-observation time of 0.1, and a full state observation likelihood model with additive Gaussian noise that has a variance of 4. The code can be run using the command matlab `sample_run.m` from the `SampleRun` folder. 

The script first defines an object called `model` that contains the parameters of the Lorenz 63 system, an array with the true hidden state, and an array with the observations to be assimilated. The script then calls the `spin_up` function to generate an initial set of representative samples from the bulk of the filtering distribution based on a spin-up phase of 2000 assimilation steps. The spin-up phase is performed with 200 samples using the perturbed observation Ensemble Kalman filter.

	L63_problem_setup
	model = spin_up(model, 200);

To assimilate the data using a default setup of the Stochastic Map filter, we define the map's structural parameters in the `options` object. These include the number of samples `M`, the pairwise distance matrix between all states `distMat`, the degree of linearity in the map `order`, the number of components before the prior-to-posterior maps revert to the identity `nonId_radius`, the spatial radius of dependence for each variable `offdiag_rad`, and the forecast sample inflation factor `rho`.

	options = struct;
	options.M            = 200;
	options.distMat      = metric_lorenz(model.d);
	options.order_all    = 2;
	options.nonId_radius = model.d;
	options.offdiag_rad  = model.d;
	options.rho          = 0.1;

Given the `model` and `options` object, the `SM` object is defined and the sequential inference is performed using the `sample_posterior` function in the `StochasticMapFilter` class. The returned `filter` variable contains the estimate for the mean, variance, covariance and quantiles of the filtering distribution over the assimilation steps.

	SM = StochasticMapFilter(model, options);
	filter = seq_assimilation(model, @SM.sample_posterior);

The script also produces a plot that compares the mean of the filtering distribution to the true hidden state for all three state variables over 1000 steps. Lastly, we also provide scripts that compare the performance of increasing non-linearity in the map and different ensemble sizes for the Lorenz 63 and Lorenz 96 systems in the `Lorenz63` and `Lorenz96` folders, respectively.
