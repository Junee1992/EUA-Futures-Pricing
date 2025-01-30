%% The density of multivariate Generalised Laplace distribution
% Kozuboski, Podgorski, Rychlik, 2013.
% x and mu are column vectors (n x 1), sigma is (n x n).
%
% Inputs:
% x : a single vector of MGL sample
% mu : a vector of means of MGL
% Sigma : a covariance matrix of MGL
% s : parameter of the standard gamma distribution
%
% dens : the density of MGL distribution

function dens = mvLaplacedensity(x, mu, Sigma, s)
    
    if ~ismatrix(x)
        x = x(:); % Ensure x is a column vector for single observations
    end

    if nargin < 2 || isempty(mu) % assume zero mean if not specified.
        mu = zeros(size(x,1),1);
    end

    if nargin < 3 || isempty(Sigma) % assume identity matrix if not specified.
        Sigma = eye(size(x, 1));
    end

    if ~ismatrix(Sigma)
        Sigma = diag(Sigma);
    end

    % Ensure Sigma is symmetric and positive-definite
    [R, p] = chol(Sigma);
    if p ~= 0
        error('Matrix Sigma is not positive-definite.');
    end

    if nargin < 4 || isempty(s) % assuming exponential (1) for W.
        s = 1;
    end

    k = size(Sigma, 1);
    Omega = inv(Sigma); % Compute the inverse of Sigma
    Q = sqrt(x'*Omega*x); C = sqrt(2+mu'*Omega*mu); % Q and C in (5).
    logdetSigma = 2 * sum(log(diag(R))); % Compute log determinant of Sigma

    dens = log(2) + mu' * Omega * x - (k/2) * log(2*pi) - log(gamma(s))- ...
        0.5 * logdetSigma + (s-k/2) * (log(Q) - log(C)) + ...
        log(besselk(s-k/2, Q*C)); % The pdf in Equation (4)
    dens = exp(dens);
end
