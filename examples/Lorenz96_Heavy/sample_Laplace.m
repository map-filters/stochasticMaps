function yL = sample_Laplace(x, theta)
% Sample elementwise from Laplace distribution by 
% generating double exponential sample (using 
% symmetry of the exponential distribution)

[N,d] = size(x);

yE = exprnd(theta,N,d);
rB = binornd(1,0.5,N,d);
yL = x(:) + yE(:).*(2*rB(:) - 1);

end
