%% Density of MGH Distribution
% dens = dghypmv(x, lambda, psi, chi, mu, sigma, gamma)
% 
% Inputs:
% x : a column vector for a sample from MGH distribution
% lambda, psi, chi : parameters of GIG distribution
% mu : a column vector of means
% sigma : a covariance matrix
% gamma : a vector of skewness parameters

function dens = dghypmv(x, lambda, psi, chi, mu, sigma, gamma)

[d,n] = size(x);
detSigma = det(sigma);
invSigma = inv(sigma);

% Calculate the Mahalanobis distance
delta = x - mu;
Q = (delta)' * invSigma * delta;

% Check if gamma is zero
if all(abs(gamma) == 0)
    symm = true;
    skewnessScaled = 0;
    skewnessNorm = 0;
else
    symm = false;
    skewnessScaled = (x - repmat(mu, n, 1))' * invSigma * gamma;
    skewnessNorm = gamma' * invSigma * gamma;
end

% Calculate the density
logTop = log(besselk(lambda - d/2, sqrt((psi + skewnessNorm) .* (chi + Q)))) + skewnessScaled;
logBottom = (d/2 - lambda) .* log(sqrt((psi + skewnessNorm) .* (chi + Q)));
logConstTop = -lambda/2 .* log(psi * chi) + (d/2) .* log(psi) + ...
    (d/2 - lambda) .* log(1 + skewnessNorm./psi);
logConstBottom = d/2 .* log(2 * pi) + log(besselk(lambda, sqrt(chi * psi))) + 0.5 .* log(detSigma);
dens = logConstTop + logTop - logConstBottom - logBottom;
dens = exp(dens);
end


