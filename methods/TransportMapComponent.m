% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef TransportMapComponent
	% Define a component of the TransportMap object that is a part
	% of the lower-triangular inverse transport maps S(x)
	% with a monotone RBF parametrization for S_{k}(x)
	%
	% Methods include: setIdComp, evaluate, inverse, Optimize

	properties

		d         % dimension of target distribution
		order     % (dx1) vector with orders for each variable in component
		OffD      % handle for non-monotone approximation object
		Diag      % handle for monotonic approximation object
        OffD_nvar % number of variables in OffD component
        Diag_nvar % number of variables in Diag component
        lambda    % L2 regularization parameter for opt
        delta     % small non-zero gradient penalty for opt
        
    end

	methods
		function TMc = TransportMapComponent(d,order,options,varargin)

			% Define TM object
			p = ImprovedInputParser;
			addRequired(p,'d');
			addRequired(p,'order');
			parse(p,d,order,varargin{:});
			TMc = passMatchedArgsToProperties(p, TMc);

            % check sizes of orders
            if length(order) ~= d
                error('Dimenion mismatch for variable orders.')
            end
			if size(order,1) ~= 1
			    error('order should be a row vector.')
			end

			% define components in TM
			TMc.OffD = NonMonotonePart(d-1, order(1:d-1));
			TMc.Diag = MonotonePart(1, order(d));
            
            % define variables
            TMc.OffD_nvar = d-1;
            TMc.Diag_nvar = 1;

            % set parameters for basis: scalingRbf, npoints_interp, kappa
            TMc.OffD.scalingRbf     = options.scalingWidths;
            TMc.Diag.scalingRbf     = options.scalingWidths;
            TMc.Diag.npoints_interp = options.npoints_interp;
            TMc.Diag.kappa          = options.kappa;
            
            % load paramters for optimization
            TMc.lambda = options.lambda;
            TMc.delta  = options.delta;
            
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function coeffs = get_coeffs(TMc)
			coeffs = [TMc.OffD.coeffs; TMc.Diag.coeffs];
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function d = ncoeff(TMc)
			d = TMc.OffD.ncoeff() + TMc.Diag.ncoeff();
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function TMc = set_Id_comp(TMc)
            TMc.OffD = TMc.OffD.set_zero_function();
			TMc.Diag = TMc.Diag.set_Id_function();
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function Sx = evaluate(TMc, X)

			% check dimension of input samples X
			if size(X,2) ~= TMc.d
				error('Map and input sample dimension mismatch')
			end
	
			% Evaluate component of transport map at samples X
            Sx = TMc.OffD.evaluate(X(:,1:end-1)) + TMc.Diag.evaluate(X(:,end));            

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function Xk = inverse(TMc, X, Z)

			% check dimension of input samples X and Z
            if size(X,2) ~= (TMc.d-1)
                error('Component and input sample dimension mismatch')
            end
            if size(Z,2) ~= 1
                error('Too many dimensions in output samples')
            end

			% evaluate OffD component at samples X(1:d-1) and invert diagonal
			SOffDx = TMc.OffD.evaluate(X);
			Xk = TMc.Diag.inverse(Z - SOffDx);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function dxS = grad_x(TMc, X)

            % check dimension of input samples X
            if size(X,2) ~= TMc.d
                error('Component and sample dimension mismatch')
            end
            
            % evaluate derivatives of off-diagonal and diagonal components
            SOff_grad_x  = TMc.OffD.grad_x(X(:,1:end-1));
            SDiag_grad_x = TMc.Diag.grad_x(X(:,end));
            
            % Evaluate derivative of k-th component of transport map
            dxS = [SOff_grad_x, SDiag_grad_x];

        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function d2xS = hess_x(TMc, X)

            % check the dimension of input samples X
            if size(X,2) ~= TMc.d
                error('Component and sample dimension mismatch')
            end
            
            % declare matrix to store hessian outputs
            N = size(X,1);
            d2xS = zeros(N, TMc.d, TMc.d);

            % save derivatives of off-diagonal and diagonal components
            d2xS(:,1:TMc.d-1,1:TMc.d-1) = TMc.OffD.hess_x(X(:,1:end-1));
            d2xS(:,TMc.d,TMc.d) = TMc.Diag.hess_x(X(:,end));
            
        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function TMc = optimize(TMc, X)

			% check the input samples
            N = size(X,1);
            if size(X,2) ~= TMc.d
                error('Check dimension of samples.')
            end
           
            % --- Setup matrices for nonlinear regression ---
           
			% evaluate monotone basis functions
            TMc.Diag = TMc.Diag.construct_basis(X(:,end));
            Psi_mon  = TMc.Diag.basis_eval(X(:,end));
            dPsi_mon = TMc.Diag.basis_grad(X(:,end));
                        
  			% Normalize monotone basis functions
			meanPsi  = sum(Psi_mon,1)/N; 
			stdPsi   =  sqrt( sum((Psi_mon-meanPsi).^2,1)/N );
			Psi_mon  = (Psi_mon - meanPsi)./stdPsi;
			dPsi_mon = dPsi_mon./stdPsi;
            
            if TMc.OffD_nvar ~= 0
                
                % evaluate off-diagonal basis functions
			    [TMc.OffD] = TMc.OffD.construct_basis(X(:,1:end-1));
                Psi_offD = TMc.OffD.basis_eval(X(:,1:end-1));
                
			    % normalize off-diagonal covariates
			    meanX = sum(Psi_offD,1)/N;  
			    stdX = sqrt( sum((Psi_offD-meanX).^2,1)/N );  
			    Psi_offD = (Psi_offD-meanX)./stdX;

                % assemble reduced QR to solve least squares problem
                [Q, R] = qr([Psi_offD; sqrt(TMc.lambda)*eye(size(Psi_offD,2))],0);
                Q1 = Q(1:N,:);
                Asqrt = Psi_mon - Q1*(Q1'*Psi_mon);
                A = (Asqrt'*Asqrt)/N;

            else
                % use only diagonal component if OffD = 0
			    A = (Psi_mon'*Psi_mon)/N;
            end

			% --- Approximate diagonal component ---
            
			% for linear S_{k}^{2} (i.e., order 0 function), use closed
			% form solution for linear monotone component
            if TMc.Diag.order == 1
                if ~isscalar(A)
			        error('Quadratic matrix should be a scalar.')
                end
			    g_mon = sqrt(1/A);

            % for nonlinear diagonal, use Newton solver for coefficients
            else
			    x0 = ones(TMc.Diag.ncoeff,1);
			    xopt = TMc.setup_and_run_optim(A, TMc.lambda, dPsi_mon, TMc.delta, x0);
			    g_mon = xopt;
            end
            
            % --- Approximate off-diagonal components ---

            if TMc.OffD_nvar ~= 0
                % compute off-diagonal coefficients
                g_off = -R\(Q1'*Psi_mon*g_mon);    
			    
			    % rescale coefficients
			    g_mon = g_mon./stdPsi(:);
			    g_off = g_off./stdX(:);
                
                % compute constant term (expectation of selected terms)
			    constTerm = -meanPsi*g_mon - meanX*g_off;
                
            else  
			    %Rescale coefficients
			    g_mon = g_mon./stdPsi(:);
                
                % compute constant term (expectation of selected terms)
			    constTerm = -meanPsi*g_mon;
                
            end
                        
            % assign coefficients to TMc.OffD and TMc.Diag
            counter = 1;
            for kk = TMc.OffD.activeVars
                order_kk = TMc.OffD.order(kk);
                TMc.OffD.coeffs{kk} = g_off(counter:counter+order_kk-1);
                counter = counter+order_kk;
            end
            TMc.OffD.constTerm = constTerm;
            TMc.Diag.coeffs = g_mon;

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function xopt = setup_and_run_optim(TMc, A, lambda, dPsi_mon, delta, x0)

			% Setup optimization
			nbasis_mon = size(A,1);
			nsamples = size(dPsi_mon,1);
			A = A + (lambda/nsamples)*eye(nbasis_mon); % add L2 regularization
			b = delta*sum(A,2);

			% run numerical optimization
			J = @(x) TMc.objectiveOptim(x, A, b, dPsi_mon, delta);
			[xopt, ~, ~, ~, ~] = projectedNewton(x0, J);
			xopt = xopt + delta;
			
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function [fx, gx, Hx]=objectiveOptim(TMc, x, A, b, dPsi, delta)
		% setup objective, gradient and hessian information for an
        % optimization problem of the form min_{a} E[0.5xAx - log(dAx + delta)]

			if size(x,2)~=1
				error('Number of coefficients should be a column vector.')
			end
			N = size(dPsi,1);

			% compute components of objective function
			Ax = A*x;
			dS = dPsi*x + delta*sum(dPsi,2);
			logdS = log(dS);

			% compute objective
			fx = .5*dot(x,Ax) - sum(logdS)/N + dot(x,b);

			% compute gradient vector
			if nargout>=2
				dPsi_dS = dPsi./dS;
				gx = Ax - sum(dPsi_dS,1)'/N + b; 
			end

			% compute hessian matrix
			if nargout>=3
				Hx = A + (dPsi_dS'*dPsi_dS)/N;
			end

        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
	end %endMethods
end %endClass