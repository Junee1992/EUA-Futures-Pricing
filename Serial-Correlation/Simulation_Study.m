%% Simulation study of Schwartz-Smith Two-Factor model with serially correlated measurement errors

% This code runs the simulation study of the two-factor model, with AR(1)
% measurement errors.
% i.e. y_t = d_t + B_t * x_t + v_t,
% v_t = M * v_{t-1} + epsilon_t,
% where M is a diagonal matrix of AR coefficients.
close all; clear all;

% Setting relevant inputs
T = [(1:5)/12]; % Maturities, say 1 - 5 months
ncontracts = size(T,2); % Number of contracts to be simulated
nobsn = 500; % Number of observations to be simulated
deltat = 1/360; % Time difference between each case.

% Model selection
ms_error = 1;
% ms_error = 1 for no independent volatilities in measurement errors,
% ms_error = 2 for inter-correlated volatilities in measurement errors.

% Model parameters (Set true values)
kappa = 2;
sigmachi = 0.1;
lambdachi = 0.01;
gamma = 1;
mu = 0.5;
sigmaxi = 0.1;
lambdaxi = 0.01;
rho = 0.8;
vsig = repelem(0.01, ncontracts); % Volatilities of measurement errors.
correl = repelem(0.8, ncontracts); % intercorrelations between measurement errors.
m = repelem(0.9, ncontracts); % AR coefficients. Set this variable as a zero vector
%if one wishes to simulate serially uncorrelated measurement errors.

par = [kappa, sigmachi, lambdachi, gamma, mu, sigmaxi, lambdaxi, rho, vsig, m];

% Choose appropriate number of parameter set.
if ms_error == 2
    par = [par, correl];
end

[Y, x, ttm] = simulateY_ar1(par, T, ncontracts, nobsn, deltat, ms_error);

% Y is the simulated price
% x is the simulates state variables
% ttm is the time to maturity

% Parameter Estimation
% Set initial values for parameters
kappa_init = 1.5;
sigmachi_init = 0.05;
lambdachi_init = 0.02;
gamma_init = 0.5;
mu_init = 1;
sigmaxi_init = 0.05;
lambdaxi_init = 0.02;
rho_init = 0.9;
vsig_init = repelem(0.1, ncontracts); % Volatilities of measurement errors.
m_init = repelem(0.6, ncontracts); % AR coefficients
correl_init = repelem(0.9, ncontracts); % Intercorrelations

% Setting lower and upper bounds for each parameter
% These boundaries can be selected arbitrarily, but with care.
lb(1:8) = [0.0001, 0.0001, 0.0001, 0.0001, -4, 0, 0, -0.9999];
lb(9:8+ncontracts) = 0.000001;
lb(9+ncontracts:8+2*ncontracts) = -0.9999;
ub(1:8) = [3, 1, 0.9999, 3, 3, 1, 0.9999, 0.9999];
ub(9:8+ncontracts) = 1;
ub(9+ncontracts:8+2*ncontracts) = 0.9999;


% Select your model
model = "price";
% Select model = "price" for modelling of logarithmic price.
% Select model = "return" for modelling of logarithmic return.

ms_error = 2;
% ms_error = 1 for no independent volatilities in measurement errors,
% ms_error = 2 for inter-correlated volatilities in measurement errors.

AR_process = "AR";
% Selected AR_process = "AR" to consider AR process in measurement errors,
% and for estimation of AR coefficients.
% Otherwise, use AR_process = "noAR" for serially uncorrelated measurement errors.

% Parameter estimation contraint to avoid parameter identification problem
% Modify boundaries accordingly for each model
if AR_process == "noAR"
    A = [-1, 0, 0, 1, 0, 0, 0, 0, repelem(0, ncontracts * 2)];
    b = 0;
    par_init = [kappa_init, sigmachi_init, lambdachi_init, gamma_init, mu_init, sigmaxi_init, lambdaxi_init, rho_init, vsig_init, correl_init];
elseif AR_process == "AR"
    lb(9+2*ncontracts:8+3*ncontracts) = -0.9999;
    ub(9+2*ncontracts:8+3*ncontracts) = 0.9999;
    A = [-1, 0, 0, 1, 0, 0, 0, 0, repelem(0, ncontracts * 3)];
    b = 0;
    par_init = [kappa_init, sigmachi_init, lambdachi_init, gamma_init, mu_init, sigmaxi_init, lambdaxi_init, rho_init, vsig_init, correl_init, m_init];
end

options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 1000, 'MaxFunEvals', 50000);
global save_att save_ytt save_ett

% Initial distribution of state variables
% E(x_{1|0}) = a0
% Cov(x_{1|0}) = P0

if model == "price"
    % As suggested in Binkowski et al (2009)
    a0 = [0 ; mu / gamma];
    P0 = [sigmachi^2 / (2*kappa), sigmachi * sigmaxi * rho / (kappa + gamma);
        sigmachi * sigmaxi * rho / (kappa + gamma), sigmaxi^2 / (2 * gamma)];
    
elseif model == "return"
    % Choose to use the first two observations for initial expectation.
    a0 = [0; Y(2,1); 0; Y(1,1)];
    P0 = [0.1, 0.1, 0, 0; 0.1, 0.1, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0];
    
end

[par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@KFS_ar, par_init, A, b, [], [], lb, ub, [], options, Y, deltat, ttm, model, a0, P0, ms_error, AR_process);

figure;
subplot(1,2,1);
plot(save_att);
hold on
plot(x);
hold off
legend("estimated \chi_t", "estimated \xi_t", "simulated \chi_t", "simulated \xi_t")

subplot(1,2,2);
plot(save_ytt(:,1));
hold on
plot(Y(:,1));
hold off
legend("estimated y", "simulated y")

for i = 1:ncontracts
    RMSE(1,i) = i;
    RMSE(2,i) = sqrt(mean(save_ett(:,i)).^2);
end
sprintf('In-sample RMSE for %d is %d.\n', RMSE)
