function y = quantiles_sorted_vector(v,p)
%~~ Copyright 1993-2010 The MathWorks, Inc. 

if ~isvector(v)
    error('should be a vector')
end
v = v(:);
n = length(v);

if isscalar(p) && (p>1)
    p = (1:p) / (p+1);  %p uniform quantiles on the ensemble
    p = p(:);
elseif isvector(p)
    p = p(:);
    if any(p < 0 | p > 1)
        error('invalid quantiles')
    end   
else
    error('case not supported')
end

% Form the vector of index values (numel(p) x 1)
r = p*n;
k = floor(r+0.5); % K gives the index for the row just before r
kp1 = k + 1;      % K+1 gives the index for the row just after r
r = r - k;        % R is the ratio between the K and K+1 rows


% Find indices that are out of the range 1 to n and cap them
k(k<1) = 1;
kp1 = min(kp1, n);

% Use simple linear interpolation for the valid percentages
y = (0.5+r).*v(kp1,:)+(0.5-r).*v(k,:);

y = reshape(y,1,length(y)); %returns always a row vector
end