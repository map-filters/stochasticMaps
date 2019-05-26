%---------- Projected Newton algorithm for strictly convex problems ----------%
function [xopt, message, rdeltaJ, norm_Pgk, it] = projectedNewton(x0, J, type, options)
%Default parameters
rtol_Jdef = 1e-6;
rtol_gdef = 1e-6;
itmaxdef = 30;
epsilon = 0.01;
switch nargin
    case 2
    rtol_J = rtol_Jdef;
    rtol_g = rtol_gdef;
    itmax = itmaxdef;
    type = 'trueHessian'; % possible types: (1) trueHessian, (2) modHessian, (3) gradient
    case 3
    rtol_J = rtol_Jdef;
    rtol_g = rtol_gdef;
    itmax = itmaxdef;    
        
    case 4 %options provided
    try
        rtol_J = options.rtol_J;
    catch
        rtol_J = rtol_Jdef;
    end
    try
        rtol_g = options.rtol_g;
    catch
        rtol_g = rtol_gdef;
    end
    try
        itmax = options.itmax;
    catch
        itmax = itmaxdef;
    end     
end
% ------- Initialize parameters -------
if size(x0,2)~=1
   error('must be a column vector') 
end
if any(x0<0)
   error('initial conditions are not feasible')
end
xk = x0; 
[Jk, gk, Hk] = J(xk);
dim = length(x0);
if ~isscalar(Jk)
    error('must be a scalar')
end
if size(gk,1)~=dim || size(gk,2)~=1
    error('wrong dimension')
end
if size(Hk,1)~=dim || size(Hk,2)~=dim
   error('wrong dimension') 
end
norm_Pg0 = norm(projectGradient(xk,gk));
tol_g = norm_Pg0*rtol_g;
norm_Pgk=norm_Pg0;
rdeltaJ=rtol_J+1;
Jold = Jk;
it = 0;
% ------- Iteration -------
while (rdeltaJ > rtol_J) && (norm_Pgk>tol_g) && (it < itmax)
    %Define search direction
    wk = norm(xk - max(0, xk - gk));
    epsk = min(epsilon, wk);
    Ik = find((xk<=epsk).*(gk>0)); %could be empty
    if ~isempty(Ik) %if Ik is not empty
        hxk = diag(Hk);
        zk = zeros(dim,1);
        zk(Ik) = hxk(Ik);
        Hk( Ik, : ) = 0;
        Hk( :, Ik ) = 0;
        Hk = Hk + diag( zk );   
    end
    switch type %Define search direction
        case 'trueHessian'
            Lk = chol(Hk,'lower'); %raises an error if Hk is not pd 
            pk =  transpose(Lk)\(Lk\gk); %search direction (== Hk\gk)

        case 'modHessian'
            itmax = 100;
            try Lk = chol(Hk,'lower'); %standrard iteration for convex problems
                pk =  transpose(Lk)\(Lk\gk);
            catch %Hessian not pd
                %Strategy 1 - flip negative eigenvalues
                if ~issymmetric(Hk)
                    Hk = (Hk + Hk')/2;
                end
                disp('hessian not pd')
                [eigVec, eigVal] = eig(Hk, 'vector');
                eigVal = abs(eigVal);
                eigVal = max(eigVal, 1e-8);
                pk =  eigVec*( (1./eigVal).*(eigVec'*gk) );   
            end
            
        case 'gradient' %revert to gradient descent
            itmax = 1000;
            pk = gk;
        otherwise
            error('method not implemented yet')
    end
    %Line search
    alphak = Armijo(xk, gk, pk, Jk, Ik, J);
    %Update
    xk = max(0, xk - alphak*pk); %iterate remains feasible
    [Jk, gk, Hk] = J(xk);  
    %Convergence criteria
    rdeltaJ = abs(Jk - Jold)/abs(Jold);
    Jold = Jk;
    norm_Pgk = norm(projectGradient(xk,gk));
    it = it + 1;
end
% ------- Exit -------
xopt = xk;
if it < itmax
    message = 'success';
else
    message = 'maxIt';
    fprintf('Reached max number of iterations during optimization \n')
end
end

function Pgk = projectGradient(xk,gk)
Pgk = gk;
Pgk( (xk==0) & (gk>=0) )=0;
end

%---------- Line search ----------%
function alpha = Armijo(xk, gk, pk, Jk, Ik, J)
%[xk] current iterate
%[gk] gradient at gk
%[pk] search direction at xk
%[Jxk] objective function at xk
%[J] @ to objective function
%[Ik] set of locally optimal coordinates on the boundary
%
%
%-------- Fixed parameters --------
itmax = 15;
sigma = 1e-4;
beta = 2;
%----------------------------------
%
it = 0; %iteration counter
Jxk_alpha_lin = Jk + 1; %enter the while loop
alpha = beta; %the first iteration inside the while loop will be for alpha = 1
while(Jk < Jxk_alpha_lin && it < itmax )
    alpha = alpha/beta;
    alpha_pk = alpha*pk;
    xk_alpha = max(0, xk - alpha_pk);
    Jxk_alpha = J(xk_alpha);
    alpha_pk(Ik) = xk(Ik) - xk_alpha(Ik);
    Jxk_alpha_lin = Jxk_alpha + sigma*dot(gk, alpha_pk);
    it = it + 1;
end
if(Jk < Jxk_alpha_lin) %max number of iterations reached
    fprintf('Line search reached max number of iterations \n')
    xkh = xk;
    xkh(Ik) = 0;
    index = find((xkh>0) & (pk>0));
    if isempty(index)
        alpha = 1;
    else
        alpha = min(xk(index)./pk(index));
    end
end
end
