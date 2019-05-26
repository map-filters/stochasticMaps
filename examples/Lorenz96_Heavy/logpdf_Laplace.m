function p = logpdf_Laplace(x,y,theta)
% Evaluate log density of Laplace (double exponential)
% distribution with parameter theta

p = -log(2*theta) - abs(y - x)/theta;

end
