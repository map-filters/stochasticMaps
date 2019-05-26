% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef StochasticMapFilter
	% Define an object from the composition of TransportMaps to
    % generate posterior (from the conditional distribution Y|X)
    % samples given samples from the prior and the likelihood model

	properties

		order      % order of transport maps
        model      % struct with the likelihood and observation data
        options    % struct with options for assimilation
		TM  	   % map from samples for approximating posterior
        NonIdComp  % list of non-identity components in map

    end

	methods 
		function SM = StochasticMapFilter(model, options, varargin)

			% define SM object
			p = ImprovedInputParser;
            addRequired(p,'model');
            addRequired(p,'options');
			parse(p, model, options, varargin{:});
			SM = passMatchedArgsToProperties(p, SM);

            % define default parameters for SM
            SM = SM.default_tmap();

            % define component orders and identify non-identity components
            [SM.order, SM.NonIdComp] = SM.setup_tmap();
            
            % define TM object
            d = length(SM.order);
            SM.TM = TransportMap(d, SM.order, SM.options);

            % set remaining components to be identity function
            IdComp = setdiff(1:SM.TM.d, SM.NonIdComp);
            for Ck=IdComp
                SM.TM.S{Ck} = SM.TM.S{Ck}.set_Id_comp();
            end
            
		end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function SM = default_tmap(SM)
        % function adds recommended options for stochastic map filter

            SM.options.scalingWidths    = 2.0;
            SM.options.lambda           = 0;
            SM.options.delta            = 1e-8;
            SM.options.npoints_interp   = max(2000,2*SM.options.M);
            SM.options.kappa            = 4;

            SM.options.data_order       = SM.options.order_all;
            SM.options.offdiag_order    = SM.options.order_all;
            SM.options.diag_order_obs   = SM.options.order_all;
            SM.options.diag_order_unobs = 1;
            SM.options.locLik           = 1;

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function [order, NonIdComp] = setup_tmap(SM)
        % define the orders for the S map components based on the specified
        % options and model parameters

            % extract parameters from model and options
            d = SM.model.d;
            distMat = SM.options.distMat;

            % set off-diagonal orders for the map
            offdiag_order = SM.options.offdiag_order;
            offdiag_rad   = SM.options.offdiag_rad;
            dist_to_order = [offdiag_order*ones(1,offdiag_rad), zeros(1,d-offdiag_rad)];

            % define cell to store orders
            order = cell(d+1,1);

            % extract ordering of variables based on distance
            [~, permutation] = sort(distMat(:,1));
            
            for i=1:d
                
                % for each variable, extract distance to observed node
                node_i = permutation(i);
                dist_i_1 = distMat(1,node_i);

                % define vector to store orders
                orders_i = zeros(1,i);

                % if dist_i_1 <= nonIdcomponents, component is identity
                if dist_i_1 <= SM.options.nonId_radius

                    % define orders for off-diagonal TM component
                    for j = 1:i-1

                        % for each previous variable, extract distance to node_i
                        node_j = permutation(j);
                        dist_i_j = distMat(node_i,node_j);

                        % compute order for variable based on distance
                        orders_i(j) = dist_to_order(dist_i_j);
                        
                    end

                    % define orders for diagonal TM component
                    if i == 1 % observed node
                        orders_i(i) = SM.options.diag_order_obs;
                    else % unobserved node
                        orders_i(i) = SM.options.diag_order_unobs;
                    end

                    % append order for dependence on data and save orders
                    if (i == 1 && SM.options.locLik == 1) || (SM.options.locLik == 0)
                        order{i+1} = [SM.options.data_order, orders_i];
                    else
                        order{i+1} = [0, orders_i];
                    end

                end
            end

            % define non-identity components of map
            NonIdComp = find(cellfun(@sum,order)~=0)';

            % set remaining components to order=[0,\dots,0,1]
            IdComp = setdiff(1:d+1, NonIdComp);
            for i=IdComp
                order{i} = [zeros(1,i-1), 1];
            end
            
        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function Xpost = assimilate_scalar_obs(SM, Xpr, data_idx, Yt)
        % assimilate scalar observation of X(:,data_idx) based on the
        % observation in Yt and by permutating the samples

            % check if Yt is a scalar
            if ~isscalar(Yt)
                error('Observation should be a scalar in SM.evaluate')
            end

            %----- Permute the ensemble -----
            [~, permutation] = sort(SM.options.distMat(:,data_idx));
            Xpr = Xpr(:,permutation);
                        
            % apply multiplicative inflation to forecast samples
            Xinfl = SM.inflate(Xpr);

            % generate samples from local likelihood
            Xi = Xinfl(:,1);
            Yi = SM.model.sampleLik(Xi);

            % compute lower map components using separability
            SM.TM = SM.TM.optimize([Yi, Xinfl], SM.NonIdComp);
            
            % generate local-likelihood samples with un-inflated X samples
            Xi  = Xpr(:,1);
            Yi  = SM.model.sampleLik(Xi);
            
            % evaluate composed map T(y,x) at [Yi, Xpr] samples
            Xpost = SM.evaluate([Yi, Xpr], Yt);
                        
            %----- Inverse permutation -----
            Xpost(:,permutation) = Xpost;
            
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function Xpost = sample_posterior(SM, Xpr, Yt)
            
            % extract number of indices for assimilation
			data_idx = SM.model.data_idx;
            n_obs = length(data_idx);
                        
			% generate analysis ensemble by sequentially assimilating data
            for i = 1:n_obs
                Xpr = SM.assimilate_scalar_obs(Xpr, data_idx(i), Yt(i));
            end
            
            % return samples in Xpost
        	Xpost = Xpr;

        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function infl_X = inflate(SM, X)
        	mean_X  = mean(X,1);
            delta_X = sqrt(1 + SM.options.rho)*(X - mean_X);
            infl_X  = delta_X + mean_X;
        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function Tx = evaluate(SM, YX, Yt)
        % evaluate the composed map  T(y,x) at the samples YX = [Y, X]

            % check if Yt is a scalar
            if ~isscalar(Yt)
                error('Observation should be a scalar in SM.evaluate')
            end

            % find total number of samples
            N = size(YX,1);

            % evaluate forward map at YX samples
            EtaYX = SM.TM.eval_map(YX);

            % set observation and evaluate inverse map
            EtaYX(:,1) = repmat(Yt, N, 1);
            YtX = SM.TM.eval_inv_map(EtaYX);

            % extract posterior samples
            Tx = YtX(:,2:end);

        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
		function dTx = gradient(SM, YX, Yt)
        % compute the gradient of the composed map T at the samples X using
        % chain rule: T(y,x) = S(y*,\cdot)^{-1} \circ S(y,x)

            % find total number of samples
            N = size(YX,1);
        
            % evaluate \nabla S(y,x) for YX inputs
            dSx = SM.TM.grad_x(YX);
            
            % evaluate \nabla S(y^{*},x) for (Y^{*},T(Y,X)) inputs
            Tyx = SM.evaluate(YX,Yt);
            dSx_ysT = SM.TM.grad_x([Yt*ones(N,1), Tyx]);
            
            % evaluate \nabla_{x} S^{-1}(x)|S using inverse function 
            % theorem as inverse of (\nabla_{x} S(y^{*},T(y,x)))
            inv_dSinvx = dSx_ysT(:,2:end,2:end);
            
            % evaluate gradient of T by composition
            dTx = zeros(N, SM.TM.d-1, SM.TM.d);
            for i=1:N
                dTx(i,:,:) = squeeze(inv_dSinvx(i,:,:))\squeeze(dSx(i,2:end,:));
            end
            
        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
	end %endMethods
end %endClass