%% Generating random vectors from Multivariate Generalised Hyperbolic distribution
% This function uses Y = mu + gamma * W + sqrt(W) * Z, Z ~ MVN(0, Sigma),
% W ~ GIG(lambda, chi, psi).
% mu and gamma are column vectors, sigma is a square matrix.
%
% Inputs:
% Lambda, psi, chi : parameters of GIG distribution
% mu : a mean vector of MGH
% Sigma : a covariance matrix of MGH
% gamma : skewness parameters
%
% Output:
% y : a sample from MGH distribution

function y = randmgh(lambda, psi, chi, mu, Sigma, gamma)

% Ensure mu is a row vector
mu = reshape(mu, [],1);

% Ensure Sigma is a matrix
Sigma = reshape(Sigma, size(mu,1), size(mu,1));

% Check if Sigma is positive definite
[~,p] = chol(Sigma);
if p > 0
    error('Matrix Sigma is not positive-definite.');
end

% Number of columns in Sigma
k = size(Sigma, 2);

w = gigrnd(lambda, psi, chi);
z = mvnrnd(zeros(1, k), Sigma, 1);

y = mu + gamma .* w + sqrt(w) .* z';

end