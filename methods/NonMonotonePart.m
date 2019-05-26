% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef NonMonotonePart
    % Defines a function c(x_{1},\dots,x_{d}) that is parametrized using
    % linear functions and RBFs and is not constrained to be monotone.
    % This class is used as part of TransportMapComponent.
    %
    % Notes: order = 0  - inactive variable
    %        order = 1  - linear
    %        order >= 2 - linear + RBF

	properties
        
        d               % number of input variables
        order           % (dx1) list of orders for each variable
        activeVars      % list of activeVars from {1,\dots,d}
  
        coeffs          % {dx1} cell of coefficients for each variable
        centers         % {dx1} cell of centers for each RBF
        widths          % {dx1} cell of widths for each RBF
        constTerm       % scalar for absolute mean of function
        scalingRbf      % scaling for RBFs

        xCache          % cache of input samples for basisCache evaluations
        basisCache      % cache of basis evaluations
        gradxCache      % cache of basis gradient evaluations
        hessxCache      % cache of basis hessian evaluations

    end

	methods
		function NMP = NonMonotonePart(d,order,varargin)
			
			% Define NM object
			p = ImprovedInputParser;
			addRequired(p,'d');
            addRequired(p,'order');
			parse(p,d,order,varargin{:});
			NMP = passMatchedArgsToProperties(p, NMP);

            % check orders
            if any(order < 0)
                error('specified orders must be greater than 0 for NMP.')
            end
            
            % define active vars (non-zero order)
            NMP.activeVars = find(order);

            % define cells to store coeffs, centers, and widths
            NMP.constTerm = 0;
            NMP.coeffs  = cell(d,1);
            NMP.centers = cell(d,1);
            NMP.widths  = cell(d,1);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function d = ncoeff(NMP)
        	d = sum(NMP.order);
        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function NMP = set_zero_function(NMP)
            NMP.constTerm = 0;
			NMP.coeffs = {zeros(NMP.ncoeff(),1)};
        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function cx = evaluate(NMP, X)
        % evaluate the linear+RBF model based on the ncoeff and the cached
        % values of the basis evaluations
            if NMP.ncoeff == 0
                cx = NMP.constTerm*ones(size(X,1),1);
            else
                %if ~isequal(X, NMP.xCache)
                %    NMP.xCache = X;
                %    NMP.basisCache = NMP.basis_eval(X);
                %end
                % check the definition of coefficients
                if size(NMP.coeffs,2) ~= 1
                    error('Coefficients should be a column cell')
                end
                %cx = NMP.basisCache*cell2mat(NMP.coeffs) + NMP.constTerm;
                cx = NMP.basis_eval(X)*cell2mat(NMP.coeffs) + NMP.constTerm;
            end
        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
        function dcx = grad_x(NMP, X)
        % evaluate gradients of the linear+RBF model based on the ncoeff 
        % and the cached values of basis evaluations
            % check if inputs are empty
            if size(X,2) == 0
                dcx = zeros(size(X,1),0);
            % if model has zero coefficients - return gradient of constant
            elseif NMP.ncoeff == 0
                dcx = zeros(size(X,1),1);
            else
                %if ~isequal(X, NMP.xCache)
                %    NMP.xCache = X;
                %    NMP.gradxCache = NMP.basis_grad(X);
                %end
                % check the definition of coefficients
                if size(NMP.coeffs,2) ~= 1
                    error('Coefficients should be a column cell')
                end
                coeff_size = reshape(cell2mat(NMP.coeffs),1,1,NMP.ncoeff);
                coeff_rep  = repmat(coeff_size,[size(X,1),size(X,2),1]);
                %dcx = sum(NMP.gradxCache.*coeff_rep,3);
                dcx = sum(NMP.basis_grad(X).*coeff_rep,3);
            end
        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function d2cx = hess_x(NMP, X)
        % evaluate hessians of the linear+RBF model based on the ncoeff and
        % the cached values of basis evaluations
            % check if inputs are empty
            if size(X,2) == 0
                d2cx = zeros(size(X,1),0,0);
            elseif NMP.ncoeff == 0
                d2cx = zeros(size(X,1),1,1);
            else
                if ~isequal(X, NMP.xCache)
                    NMP.xCache = X;
                    NMP.hessxCache = NMP.basis_hess(X);
                end
                % check the definition of coefficients
                if size(NMP.coeffs,2) ~= 1
                    error('Coefficients should be a column cell')
                end
                coeff_size = reshape(cell2mat(NMP.coeffs),1,1,1,NMP.ncoeff);
                coeff_rep  = repmat(coeff_size,[size(X,1),size(X,2),size(X,3),1]);
                d2cx = sum(NMP.hessxCache.*coeff_rep,4);
            end
        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function basisEval = basis_eval(NMP, X)
            
             % find total number of samples
             N = size(X,1);

             % declare matrix to store basis evaluations
             basisEval = zeros(N, NMP.ncoeff); 

             % declare counter to store basis functions
             counter = 1;
             
             % evaluate basis functions for each active dimension
             for ii = NMP.activeVars
                
                % extract order and samples for each dimension
                order_ii = NMP.order(ii);
                samples_ii = X(:,ii);
                
                % linear basis evaluations
                if order_ii == 1
                    basisEval(:,counter) = samples_ii;
                    counter = counter + 1;
                % linear + RBF basis evaluations
                else
                    rbfEval = NMP.eval_RBF(samples_ii, NMP.centers{ii}, NMP.widths{ii});
                    basisEval(:,counter:counter+order_ii-1) = [samples_ii, rbfEval];
                    counter = counter + order_ii;
                end
                
             end
             
             % check total number of basis functions
             if (counter-1) ~= NMP.ncoeff
                error('NMP: check total number of basis functions.')
             end

        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
        function gradxEval = basis_grad(NMP, X)

             % find total number of samples
             N = size(X,1);

             % declare matrix to store gradient of basis evaluations
             gradxEval = zeros(N, NMP.d, NMP.ncoeff);
             
             % declare counter to store basis functions
             counter = 1;
             
             % evaluate basis functions for each active dimension
             for ii = NMP.activeVars
                
                % extract order and samples for each dimension
                order_ii = NMP.order(ii);
                samples_ii = X(:,ii);
                
                % linear basis evaluations
                if order_ii == 1
                    gradxEval(:,ii,counter) = ones(N,1);
                    counter = counter + 1;
                % linear + RBF basis evaluations
                else
                    rbfGrad = NMP.grad_x_RBF(samples_ii, NMP.centers{ii}, NMP.widths{ii});
                    gradxEval(:,ii,counter:counter+order_ii-1) = [ones(N,1), rbfGrad];
                    counter = counter + order_ii;
                end
                
             end

        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
        function hessxEval = basis_hess(NMP, X)
        % Note: hessian is diagonal because of seperability of expansion
        
            % find total numeber of samples
            N = size(X,1);
            
            % declare matrix to store gradient of basis evaluations
            hessxEval = zeros(N, NMP.d, NMP.d, NMP.ncoeff);
            
            % declare counter to store basis functions
            counter = 1;

            % evaluate basis functions for each active dimension
            for ii = NMP.activeVars
                
                % extract order and samples for each dimension
                order_ii = NMP.order(ii);
                samples_ii = X(:,ii);

                % linear basis evaluations
                if order_ii == 1
                    hessxEval(:,ii,ii,counter) = zeros(N,1);
                    counter = counter + 1;
                else
                    rbfHess = NMP.hess_x_RBF(samples_ii, NMP.centers{ii}, NMP.widths{ii});
                    hessxEval(:,ii,ii,counter:counter+order_ii-1) = [zeros(N,1), rbfHess];
                    counter = counter + order_ii;
                end
                
            end
        
        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
        function r = eval_RBF(NMP, X, centers, widths)
        % evaluate 1D RBFs for all centers and widths in vectors
            r = exp(-((X - centers)./widths).^2/2)./(widths*sqrt(2*pi));
        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
        function dr = grad_x_RBF(NMP, X, centers, widths)
        % evaluate grad of 1D RBFs for all centers and widths in vectors
            r  = exp(-((X - centers)./widths).^2/2)./(widths*sqrt(2*pi));
            dr = r*(-1).*((X - centers)./widths.^2);
        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
        function d2r = hess_x_RBF(NMP, X, centers, widths)
        % evaluate hess of 1D RBFs for all centers and widths in vectors
            r = exp(-((X - centers)./widths).^2/2)./(widths*sqrt(2*pi));
            d2r = r.*(((X - centers)./widths.^2).^2 - 1./widths.^2);
        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
        function NMP = construct_basis(NMP, X)
        % set the centers and widths based on the specified order for each
        % variable dimension.

            % check dimension of inputs
            if size(X,2) ~= NMP.d
                error('Dimension mismatch')
            end

            % define the RBF widths and centers for each active dimension
            for ii = NMP.activeVars

                % compute the widths and centers
                order_ii = NMP.order(ii);
                if order_ii == 1
                    centers_ii=[]; widths_ii=[];
                else
                    [centers_ii, widths_ii] = NMP.centers_and_widths(X(:,ii), order_ii-1);
                end
                
                % save widths and centers
                NMP.centers{ii} = centers_ii;
                NMP.widths{ii}  = widths_ii;
                
            end

        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function [centers, widths] = centers_and_widths(NMP, X, ncoeff)
		% compute the centers and widths of the RBF basis specified using
        % ncoeff functions using the samples X
        % 
        % Note: input X can be unsorted (sort before computing centers)
            
            if size(X,2) ~= 1
				error('Samples should be a column vector.')
            end
            if ~isscalar(NMP.scalingRbf) || NMP.scalingRbf<= 0
				error('scalingRBF be a positive scalar.')
            end
            if ncoeff < 1
				error('centers_and_widths should input at least one basis.')
            end
            
            % sort samples in X
            X = sort(X,1);

			% define the center and widths based on the quantiles
            if ncoeff == 1
			    qq = quantiles_sorted_vector(X,[.25, .5, .75]);   
			    centers = qq(2);
			    widths  = (qq(3)-qq(1))/2;
            else
                % define the centers of the RBFs
            	centers = quantiles_sorted_vector(X, ncoeff);
				% define the widths of the RBFs
				if ncoeff == 2
					widths = (centers(2)-centers(1))*ones(1,2);
				else
					widths = zeros(1,ncoeff);
					widths(2:end-1) = (centers(3:end)-centers(1:end-2))/2;
					widths(1)       = centers(2)-centers(1);
					widths(end)     = centers(end)-centers(end-1);
				end
            end
            
            % scale the widths by NMP.scalingRBF
            widths = NMP.scalingRbf*widths;

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
    end %endMethods
end %endClass