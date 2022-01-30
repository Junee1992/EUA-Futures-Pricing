%% Simulation for futures prices with serially correlated measurement errors
% The simulation is done under the setting of Schwartz-Smith Two-Factor
% model, with two state variables following O-U processes.
% Also, we assume the measurement errors are serially correlated with AR(1)
% process.

% Mathematical expressions for the linear state-space model
% x_t = C + G * x_{t-1}+ + w_t, w_t ~ N(0, W)
% y_t = d_t + B_t * x_{t-1} + v_t + M * v_{t-1}, v_t ~ N(0, V)

% For simplicity, we assume 360 days in a year, 30 days in a month.
% For model 1, we assume V is a diagonal matrix. i.e. No intercorrelation between
% contracts with different maturities.
% For model 2, we assume V is not a full covariance matrix. i.e.
% Intercorrelations exist between measurement errors of different contracts.

function [simY, x, ttm, v] = simulateY_ar1(par, T, ncontracts, nobsn, deltat, model)
% par is a vector of relevant parameters of the model.
% T is a vector of maturity time in years. For example, 1/12 = 1 month.
% ncontracts is the number of contracts, and nobsn is the number of
% observations.
% deltat is the time difference between each observation. 
%i.e. deltat = t_i - t_{i-1}.

kappa = par(1);
sigmachi = par(2);
lambdachi = par(3);
gamma = par(4);
mu = par(5);
sigmaxi = par(6);
lambdaxi = par(7);
rho = par(8);
vsig = par(9:8+ncontracts); % volatilities of measurement errors.
m1 = par(9+ncontracts:8+2*ncontracts); % AR coefficients
M = diag(m1);

% Covariance matrix of measurement error
if model == 1
    V = diag(vsig.^2);
elseif model == 2
    correl = par(9+2*ncontracts:8+3*ncontracts);
   
    % Manually creating the correlation matrix of measurement errors
    CorMat = diag(repelem(1, ncontracts));
    for i = 1:ncontracts
        for j = 1:ncontracts
            if i == j
                CorMat(i,j) = 1;
            else
                CorMat(i,j) = correl(i) * correl(j);
            end
        end
    end
    D = diag(vsig);
    V = D * CorMat * D;
end

% Time to maturity
ttm = zeros(nobsn, ncontracts);
yeardays = 360;
monthdays = 30;

for j = 1:ncontracts
    for i = 1:nobsn
        if j == 1
            ttm(i,j) = T(j) - i/yeardays;
            if ttm(i,j) < 0
                ttm(i,j) = ttm((i - monthdays * 12 * T(j)), j);
            end
        else
            ttm(i,j) = T(j) - i/yeardays;
            if i > monthdays
                ttm(i,j) = ttm(i - monthdays, j);
            end
        end
    end
end

% State Equation (x_t)
C = [0 ; (mu/gamma)*(1-exp(-gamma * deltat))];
G = [exp(-kappa * deltat), 0; 0, exp(-gamma * deltat)];
W = [(1-exp(-2*kappa*deltat))/(2*kappa)*sigmachi^2, (1-exp(-(kappa + gamma)*deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho);...
    (1-exp(-(kappa + gamma)*deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho), (1-exp(-2*gamma*deltat))/(2*gamma)*sigmaxi^2];

x = zeros(nobsn, size(C,1));
x(1,:) = [0 ; mu / gamma]'; % Initial value of x.

% Measurement Equation (y_t)
d1 = (1 - exp(-2 * kappa * ttm)) * sigmachi^2 / (2 * kappa);
d2 = (1 - exp(-2 * gamma * ttm)) * sigmaxi^2 / (2 * gamma);
d3 = (1 - exp(- (kappa + gamma) * ttm)) * 2 * sigmachi * sigmaxi * rho / (kappa + gamma);
d = (mu - lambdaxi) / gamma * (1 - exp(-gamma * ttm)) - (lambdachi/kappa) * (1 - exp(-kappa * ttm)) + (1/2) * (d1 + d2 + d3);
d = d';

B = zeros(ncontracts, size(C,1), nobsn);
for i = 1:nobsn
    B1(i,:) = exp(-kappa * ttm(i,:));
    B2(i,:) = exp(-gamma * ttm(i,:));
    B(:,:,i) = [B1(i,:); B2(i,:)]';
end

% Distribution of the state error.
normmu = repelem(0, size(x, 2));
normsigma = [1, rho; rho, 1];
burn_in = 500;

% Initial values of state variables.
chi(1) = x(1,1);
xi(1) = x(1,2);

% Simulating state variables
for i = 2:nobsn + burn_in + 1
    z(i-1, :) = mvnrnd(normmu, normsigma);
    chi(i) = chi(i-1) * exp(-kappa * deltat) + sigmachi * sqrt((1 - exp(-2 * kappa * deltat)) / (2 * kappa)) * z(i-1, 1);
    xi(i) = xi(i-1) * exp(-gamma * deltat) + mu / gamma * (1 - exp(-gamma * deltat)) + sigmaxi * sqrt((1 - exp(-2 * gamma * deltat)) / (2 * gamma)) * z(i-1, 2);
end

x = [chi;xi]';
x = x(burn_in + 2:end, :); % Removing burn-in period.

% Simulating measurable variables
vnormmu = repelem(0, ncontracts);
vnormsigma = V;
y = zeros(nobsn, ncontracts);
v = mvnrnd(vnormmu, vnormsigma, nobsn);

for i = 1:nobsn
    
    if i == 1
        y(i,:) = d(:,i) + B(:,:,i) * x(i,:)' + v(i,:)';
    else
        v(i,:) = M * v(i-1,:)' + v(i,:)';
        y(i,:) = d(:,i) + B(:,:,i) * x(i,:)' + v(i,:)';
    end
   
end

simY = y;

% ------------------------ End of simulation ----------------------------%
