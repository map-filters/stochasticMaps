% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef MonotonePart
    % Define a univariate nonlinear and monotone function that is 
    % parametrized using a linear term and RBFs.
    %
    % Notes: Function doesn't implement constant term
    %        order = 1 - linear
    %        order > 1 - (order+1) RBFs

	properties

        d           % number of input variables
        order       % (dx1) list of orders for each variable
        
        coeffs      % (ncoeff x 1) array of coefficients for each basis
        centers     % (ncoeff x 1) array of RBF centers
        widths      % (ncoeff x 1) array of RBF widths
        
        xCache      % cache of input samples where basisCachce is evaluated
        basisCache  % cache of basis evaluations
        gradxCache  % cache of grad_x of basis functions
        hessxCache  % cache of hess_x of basis functions
        
        scalingRbf       % scaling for RBFs
        kappa            % domain scaling for interpolating map
        npoints_interp   % number of points used for interpolation
  
	end

	methods
		function MP = MonotonePart(d,order, varargin)
            
            % check order (currently implemented for univariate function)
            if d ~= 1
                error('MonotonePart should be one-dimensional.')
            end
            
            % Define MP object
            p = ImprovedInputParser;
            addRequired(p,'d');
            addRequired(p,'order');
            parse(p,d,order,varargin{:});
            MP = passMatchedArgsToProperties(p, MP);
            
            % check order
            if any(order <= 0)
                error('specified order must be greater than 0 for MP.')
            end
            
		end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function d = ncoeff(MP)
            % 1 basis function (linear term) if MP.order == 1
            if MP.order == 1
                d = 1;
            % p+2 basis functions (RBFs + 2 erf functions) if MP.order > 1
            else
                d = (MP.order - 1) + 2;
            end
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function MP = set_Id_function(MP)
            % set coeff, widths, centers, to identity map
            MP.coeffs = 1;
            MP.widths = [];
            MP.centers = [];
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function cx = evaluate(MP, X)
        % evaluate the mon RBF model based on the cached basis evaluations
            %if ~isequal(X, MP.xCache)
            %    MP.xCache = X;
            %    MP.basisCache = MP.basis_eval(X);
            %end
            %cx = MP.basisCache*MP.coeffs;
            cx = MP.basis_eval(X)*MP.coeffs;
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function dcx = grad_x(MP, X)
        % evaluate the grad of monotone RBFs based on the cached basis
            %if ~isequal(X, MP.xCache)
            %    MP.xCache = X;
            %    MP.gradxCache = MP.basis_grad(X);
            %end 
            %dcx = MP.gradxCache*MP.coeffs;
            dcx = MP.basis_grad(X)*MP.coeffs;
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function basisEval = basis_eval(MP, X)
            % evaluate only linear term for order = 1
            if MP.order == 1
                basisEval = X;
            % evaluate monotone RBF terms for order > 1
            else
                basisEval = MP.eval_monotone_rbfs(X, MP.centers, MP.widths);
            end            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function gradxEval = basis_grad(MP, X)
            % evaluate only constant terms for order = 1
            if MP.order == 1
                gradxEval = ones(size(X,1),1);
            % evaluate gradient of monotone RBF terms for order > 1
            else
                gradxEval = MP.grad_x_monotone_rbfs(X, MP.centers, MP.widths);
            end
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function hessxEval = basis_hess(MP, X)
            % evaluate only constant terms for order = 1
            if MP.order == 1
                hessxEval = zeros(size(X,1),1);
            % evaluate hessian of monotone RBF terms for order > 1
            else
                hessxEval = MP.hess_x_monotone_rbfs(X, MP.centers, MP.widths);
            end
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function Xk = inverse(MP, Z)

            if MP.order == 1
                Xk = Z/MP.coeffs(1);
            else
                
                % define the domain for interpolation
                lbound = MP.centers(1) - MP.kappa*MP.widths(1);
                ubound = MP.centers(end) + MP.kappa*MP.widths(end);
                xx_interp = linspace(lbound, ubound, MP.npoints_interp)';

                % evaluate the monotone function at xx_interp
                Psi_mon = MP.eval_monotone_rbfs(xx_interp, MP.centers, MP.widths);
                yy_interp = Psi_mon*MP.coeffs;
                
                % invert the function at the points in Z
                Xk = MP.invert1D_map_interp(Z, xx_interp, yy_interp);

            end

        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
        function interp_values = invert1D_map_interp(MP, points, xx, yy)
        % given ORDERED evaluation pairs [xx],[yy] from a monotone function
        % y=f(x), invert f(x) at [points] using linear interpolation

            if size(points,2)~=1  || size(xx,2)~=1 || size(yy,2)~=1
                error('should be column vectors')
            end
            if length(xx)~=length(yy)
                error('should have the same length')
            end
            n = length(yy);

            %------ Bracketing values ---------------
            [val, indPlus] = max( points<=yy', [], 2);

            indMin = indPlus - 1;

            % Fix samples above the upper range
            out_sup = find(val == 0);
            indMin(out_sup) = n-1;
            indPlus(out_sup) = n;

            % Fix samples below the lower range
            out_inf = find(indMin < 1);
            indMin(out_inf)=1;
            indPlus(out_inf)=2;

            %------ Linear interpolation ---------------
            delta = (points - yy(indMin))./( yy(indPlus) - yy(indMin) );
            interp_values = (1-delta).*xx(indMin) + delta.*xx(indPlus);

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function MP = construct_basis(MP, X)

            % check the size of samples
            if size(X,2) ~= 1
                error('Samples should be a column vector.')
            end
            % check value for MP.order
            if MP.order < 1
                error('Monotone part order is less than 1 - not allowed.')
            end
            
            % if order == 1 - set basis to be a linear function
            % if order > 1 - sort samples and use quantiles to set centers
            if MP.order == 1
                MP.centers = []; MP.widths = [];
            else
                [MP.centers, MP.widths] = MP.centers_and_widths(X, MP.ncoeff);
            end
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function [centers, widths] = centers_and_widths(MP, X, ncoeff)
        % compute the centers and widths of the RBF basis specified using
        % ncoeff functions using the samples X
        % 
        % Note: input X can be unsorted (sort before computing centers)
 
            if size(X,2) ~= 1
				error('Samples should be a column vector.')
            end
            if ~isscalar(MP.scalingRbf) || MP.scalingRbf<= 0
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
            
            % scale the widths by MP.scalingRBF
            widths = MP.scalingRbf*widths;

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function f = eval_monotone_rbfs(MP, X, centers, widths)

            % check that X is a column vector
            if size(X,2)~=1
                error('Samples for MP should be a column vector.')
            end
            N = length(X);

            % check that widths and centers are a row vector
            if size(widths,1)~=1 || size(centers,1)~=1
                error('Widths and centers should be a row vector in MP.')
            end
            % check that the number of centers and widths match
            if length(centers)~=length(widths)
                error('Number of centers and widths should mathch in MP.')
            end
            % check that at least 3 centers are specified
            if length(centers)<3
                error('Basis should consists of at least 3 elements in MP')
            end
            
            % declare matrix to store basis values
            f = zeros(N, MP.ncoeff);
            
            % evaluate centered inputs
            deltaX = (X - centers)./widths/sqrt(2);
            
            % evaluate basis functions
            f(:,1) = .5*( sqrt(2)*widths(1)*deltaX(:,1).*(1-erf(deltaX(:,1))) ...
                - widths(1)*sqrt(2/pi)*exp(-(deltaX(:,1)).^2) );
            f(:,2:end-1) = .5*(1+erf(deltaX(:,2:end-1)));
            f(:,end) = .5*( sqrt(2)*widths(end)*deltaX(:,end).*(1+erf(deltaX(:,end))) ...
                + widths(end)*sqrt(2/pi)*exp(-(deltaX(:,end)).^2) );
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function df = grad_x_monotone_rbfs(MP, X, centers, widths)

            % check that X is a column vector
            if size(X,2)~=1
                error('Samples for MP should be a column vector.')
            end
            N = length(X);

            % check that widths and centers are a row vector
            if size(widths,1)~=1 || size(centers,1)~=1
                error('Widths and centers should be a row vector in MP.')
            end
            % check that the number of centers and widths match
            if length(centers)~=length(widths)
                error('Number of centers and widths should mathch in MP.')
            end
            % check that at least 3 centers are specified
            if length(centers)<3
                error('Basis should consists of at least 3 elements in MP')
            end
            
            % declare matrix to store basis values
            df = zeros(N, MP.ncoeff);
            
            % evaluate centered inputs
            deltaX = (X - centers)./widths/sqrt(2);            
            
            % evaluate gradient of basis functions
            df(:,1)=.5*(1-erf(deltaX(:,1)));
            df(:,2:end-1) = exp( -deltaX(:,2:end-1).^2 )./(widths(2:end-1)*sqrt(2*pi));
            df(:,end)=.5*(1+erf(deltaX(:,end)));
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function d2f = hess_x_monotone_rbfs(MP, X, centers, widths)

            % check that X is a column vector
            if size(X,2)~=1
                error('Samples for MP should be a column vector.')
            end
            N = length(X);

            % check that widths and centers are a row vector
            if size(widths,1)~=1 || size(centers,1)~=1
                error('Widths and centers should be a row vector in MP.')
            end
            % check that the number of centers and widths match
            if length(centers)~=length(widths)
                error('Number of centers and widths should mathch in MP.')
            end
            % check that at least 3 centers are specified
            if length(centers)<3
                error('Basis should consists of at least 3 elements in MP')
            end
            
            % declare matrix to store basis values
            d2f = zeros(N, MP.ncoeff);
            
            % evaluate centered inputs
            deltaX = (X - centers)./widths/sqrt(2);            
            
            % evaluate gradient of 2:end-1 basis functions
            df_mid = exp( -deltaX(:,2:end-1).^2 )./(widths(2:end-1)*sqrt(2*pi));
            
            % evaluate hessian of basis functions
            d2f(:,1)= -exp( -deltaX(:,1).^2 )./(widths(1)*sqrt(2*pi));
            d2f(:,2:end-1) = -df_mid.*deltaX(:,2:end-1)*sqrt(2)./widths(2:end-1);
            d2f(:,end)= exp( -deltaX(:,end).^2 )./(widths(end)*sqrt(2*pi));

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
    end %endMethods
end %endClass