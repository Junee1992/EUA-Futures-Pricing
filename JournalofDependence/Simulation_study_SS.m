%% Simulation study of Schwartz-Smith Two-Factor model with serially correlated measurement errors

% This code runs the simulation study of the two-factor model, with AR(1)
% measurement errors.
% i.e. y_t = d_t + B_t * x_t + v_t,
% v_t = M * v_{t-1} + epsilon_t,
% where M is a diagonal matrix of AR coefficients.
close all; clear all;

% Setting relevant inputs
T = [(1:5)]; % Maturities in years, say 1 - 5 months
ncontracts = size(T,2); % Number of contracts to be simulated
nobsn = 1000;
% For-loop used for testing simulation studies with different sample sizes in the same run.
deltat = 1/260; % Time difference between each case.

% Model selection
correlation = 0;
% correlation = 0 for independent volatilities in measurement errors,
% correlation = 1 for inter-correlated volatilities in measurement errors.

LT = "OU"; % = "GBM" if the long-term factor follows GBM, "OU" if mean-reverting.

% Model parameters (Set true values)
kappa = 0.5;
sigmachi = 0.2;
lambdachi = 0.05;
if LT == "GBM"
    gamma = [];
elseif LT == "OU"
    gamma = 0.3;
end
mu = 0;
sigmaxi = 0.2;
lambdaxi = 0.05;
rho = 0.3;
vsig = repelem(0.01, ncontracts); % Volatilities of measurement errors.
correl = repelem(0.8, ncontracts); % intercorrelations between measurement errors.
n_lags = 0; % AR order (p)
% m = repelem(0.4, n_lags * ncontracts); % AR coefficients.

par = [kappa, sigmachi, lambdachi, gamma, mu, sigmaxi, lambdaxi, rho, vsig];

% Choose appropriate number of parameter set.
if correlation == 1
    par = [par, correl];
else
    par = par;
end

if n_lags >= 1
    serial = "yes";
    par = [par, m];
else
    serial = "no";
    par = par;
end

par_names = define_parameters(LT, ncontracts, correlation, n_lags)

%  for i = 1:size(nobsn_v,2)
rng(100);
[y, x, ttm, v] = simulateY(par, par_names, T, ncontracts, nobsn, deltat, LT, correlation, serial);
[nobsn, ncontracts] = size(y);

% Define the parameters
LT = "OU"; % GBM or OU.
n_season = 0;
correlation = 0; % 0 for diagonal matrix, 1 for full matrix.
max_lags = 0; % number of lags to be considered for serial correlation of measurement errors
n_lag = 0; % Assume no AR in the first iteration.

detrend_price = "no"; % "yes" for deseasonalisation, otherwise "no"
par_names = define_parameters(LT, ncontracts, correlation, n_lag)';
n_par = size(par_names, 1);

% Log-likelihood
% log_L = kf_v1(par_init, par_names, y, deltat, ttm, LT, correlation, serial, test)
% log_L = kf_v2(par_init, par_names, y, deltat, ttm, LT, correlation, serial)

global save_ett save_ytt save_att save_vt
[par_optimised, log_L, par_init, trend, season] = param_estim(y, ttm, deltat, detrend_price, n_par, ncontracts, par_names, LT, correlation, max_lags);


