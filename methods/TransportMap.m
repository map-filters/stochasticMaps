% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef TransportMap
    % Define a transport map S of dimension d that pushes forward
    % a target density specified by samples to a standard normal

    properties

        d       % number of total dimensions
        order   % {d x 1} cell of orders for each map component
        S       % {d x 1} cell of map components

    end

    methods  
        function TM = TransportMap(d, order, options, varargin)

            % Define TM object
            p = ImprovedInputParser;
            addRequired(p,'d');
            addRequired(p,'order');
            parse(p,d,order,varargin{:});
            TM = passMatchedArgsToProperties(p, TM);

            % define TM components
            TM.S = cell(TM.d,1);
            for k=1:TM.d
                TM.S{k} = TransportMapComponent(k,order{k}, options);
            end

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function TM = optimize(TM, X, nonIdComp)
        % optimize: Function estimates parameters of transport
        % map independently for each component from samples of \pi
        % using the maximum likelihood method.
        %
        % Inputs: TM - transport map object
        %         X - N x d data matrix
        %         nonIdComp - list of non-identity components

            % estimate each component using seperability of objective
            for Ck=nonIdComp
                TM.S{Ck} = TM.S{Ck}.optimize(X(:,1:Ck));
            end

            % set remaining components to be identity function
            %IdComp = setdiff(1:TM.d, nonIdComp);
            %for Ck=IdComp
            %    TM.S{Ck} = TM.S{Ck}.set_Id_comp();
            %end

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function Z = eval_map(TM, X)
        % eval_map: apply map to samples X~\pi. If S(x) is true map,
        % function maps samples from \pi to samples from \rho
        % reference distribution (i.e., standard Gaussian)
        %
        % Inputs: X - (N x d) sample matrix

            % generate matrix to store output samples
            Z = zeros(size(X));

            % evaluate Sk and save samples in Z(:,k)
            for k=1:TM.d
                Z(:,k) = TM.S{k}.evaluate(X(:,1:k));
            end

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function X = eval_inv_map(TM, Z)
        % eval_inv_map: apply map to samples Z~N(0,I). If S(x) is true map,
        % function maps samples from N(0,1) to samples from \pi (target)

            % generate matrix to store output samples
            X = zeros(size(Z));

            % invert each map component by running 1D root-solves
            for k=1:TM.d
                X(:,k) = TM.S{k}.inverse(X(:,1:k-1), Z(:,k));
            end

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function dSx = grad_x(TM,X)
        % grad_x: compute the gradients of each component function with
        % respect to the N samples in X. The output is size (N x d x d)
            
            % generate matrix to store gradients
            N = size(X,1);
            dSx = zeros(N, TM.d, TM.d);
            
            % compute the gradient for each component seperately
            for k=1:TM.d
                dSx(:,k,1:k) = TM.S{k}.grad_x(X(:,1:k));
            end
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
    end %endMethods
end %endClass