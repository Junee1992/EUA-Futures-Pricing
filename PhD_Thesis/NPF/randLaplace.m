%% Generating random vectors from Multivariate Laplace distribution
% [x, w] = randLaplace(mu, Sigma, s, n)
% This function uses Y = mu * W + sqrt(W) * Z, Z ~ MVN(0, Sigma), W ~ exp(1).
% mu is a column vector, sigma is a square matrix.
%
% Inputs:
% mu = a vector of means
% Sigma = covariance matrix
% s = parameter of standard gamma distribution in W
% n = number of samples to be generated
%
% Outputs:
% x = Samples from MGL distribution
% w = Standard gamma samples used to sample x.

function [x,w] = randLaplace(mu, Sigma, s, n)
    
    % Ensure mu is a row vector
    mu = reshape(mu, 1, []);
    
    % If Sigma is missing, create it as an identity matrix of size mu
    if nargin < 3
        Sigma = eye(length(mu));
    end

    
    % Ensure Sigma is a matrix
    Sigma = reshape(Sigma, size(mu,2), size(mu,2));
    
    % Check if Sigma is positive definite
    [~,p] = chol(Sigma);
    if p > 0
        error('Matrix Sigma is not positive-definite.');
    end
    
    % Number of columns in Sigma
    k = size(Sigma, 2);
    
    % Repeat mu if n is greater than the number of rows in mu
    if n > size(mu, 1)
        mu = repmat(mu, n, 1);
    end
    
    % Generate exponential and normal distributions
    w = gamrnd(s, 1, n, 1);
    z = mvnrnd(zeros(1, k), Sigma, n);
    
    % Calculate x
    x = mu .* w + sqrt(w) .* z;
end
