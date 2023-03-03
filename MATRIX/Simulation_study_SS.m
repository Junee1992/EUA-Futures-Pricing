%% Simulation study of Schwartz-Smith Two-Factor model with serially correlated measurement errors

% This code runs the simulation study of the two-factor model, with AR(p)
% measurement errors.
% i.e. y_t = d_t + B_t * x_t + v_t,
% v_t = phi_1 * v_{t-1} + phi_2 * v_{t-2} + ... + phi_p * v_{t-p} + epsilon_t,
% where phi_p are diagonal matrices of AR coefficients at order p.
close all; clear all;

%% Simulation Study 1

% Setting relevant inputs
T = [(1:5)]; % Maturities in years, say 1 - 5 months
ncontracts = size(T,2); % Number of contracts to be simulated
nobsn = 1000;
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
mu = 1;
sigmaxi = 0.2;
lambdaxi = 0.05;
rho = 0.3;
vsig = repelem(0.01, ncontracts); % Volatilities of measurement errors.
correl = repelem(0.8, ncontracts); % intercorrelations between measurement errors.
n_lags = 1; % AR order (p)
m = repelem(0.4, n_lags * ncontracts); % AR coefficients.

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

par_names = define_parameters(LT, ncontracts, correlation, n_lags)';
    rng(100); % for N=5, AR(1)
    
[y, x, ttm, v] = simulateY(par, par_names, T, ncontracts, nobsn, deltat, LT, correlation, serial);
[nobsn, ncontracts] = size(y);

% Define the parameters for estimation
LT = "OU"; % GBM or OU.
n_season = 0;
correlation = 0; % 0 for diagonal matrix, 1 for full matrix.
max_lags = 1; % number of lags to be considered for serial correlation of measurement errors
detrend_price = "no"; % "yes" for deseasonalisation, otherwise "no"

global save_ett save_ytt save_att save_vt

% For two-step approach: we start by asssuming no AR in the measurement
% error
par_names_model1 = define_parameters(LT, ncontracts, correlation, 0)';
n_par = length(par_names_model1);
[par_model2, log_L_model1, par_init, trend, season, att, ytt, ett] = param_estim(y, ttm, deltat, detrend_price, n_par, ncontracts, par_names_model1, LT, correlation, max_lags);
par_model1 = par_model2(1:8+ncontracts);
AIC_model1 = 2*log_L_model1+ 2*(length(par_names_model1));

serial = "yes";
par_names_model2 = define_parameters(LT, ncontracts, correlation, 1)';
log_L_model2 = kf_v2_arp(par_model2, par_names_model2, y, deltat, ttm, LT, correlation, serial);
AIC_model2 = 2*log_L_model2 + 2*(length(par_names_model2));

n_lag = max_lags;
par_names_model3 = define_parameters(LT, ncontracts, correlation, n_lag)';
n_par = length(par_names_model3);
par_init = par_model2;
[par_model3, log_L_model3, par_init, trend, season, att, ytt, ett] = param_estim_arp(y, ttm, deltat, detrend_price, n_par, par_names_model3, LT, correlation, par_init);
AIC_model3 = 2*log_L_model3 + 2*(length(par_names_model3));

model1_result = set_parameters(LT, ncontracts, par_names_model1, par_model1, correlation, "no");
model2_result = set_parameters(LT, ncontracts, par_names_model2, par_model2, correlation, "yes");
model3_result = set_parameters(LT, ncontracts, par_names_model3, par_model3, correlation, "yes");

%% Simulation Study 2

T = [(1:5)]; % Maturities in years, say 1 - 5 months
ncontracts = size(T,2); % Number of contracts to be simulated
nobsn = 1000;
deltat = 1/260; % Time difference between each case.

% Model selection
correlation = 1;
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
mu = 1;
sigmaxi = 0.2;
lambdaxi = 0.05;
rho = 0.3;
vsig = repelem(0.01, ncontracts); % Volatilities of measurement errors.
correl = repelem(0.8, ncontracts); % intercorrelations between measurement errors.
n_lags = 2; % AR order (p)
m = repelem(0.4, n_lags * ncontracts); % AR coefficients.

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

% par_names = define_parameters(LT, ncontracts, correlation, n_lags);

% for i = 1:size(nobsn_v,2)
par_names = define_parameters(LT, ncontracts, correlation, n_lags)';
    rng(106); %for N=5, AR(2), full covariance
    
[y, x, ttm, v] = simulateY(par, par_names, T, ncontracts, nobsn, deltat, LT, correlation, serial);
[nobsn, ncontracts] = size(y);

% Define the parameters
LT = "OU"; % GBM or OU.
n_season = 0;
correlation = 1; % 0 for diagonal matrix, 1 for full matrix.
max_lags = 2; % number of lags to be considered for serial correlation of measurement errors
n_lag = 0; % Assume no AR in the first iteration.

detrend_price = "no"; % "yes" for deseasonalisation, otherwise "no"
par_names = define_parameters(LT, ncontracts, correlation, n_lag)';
n_par = size(par_names, 1);

global save_ett save_ytt save_att save_vt

LT = "OU"; % GBM or OU.
n_season = 0;
correlation = 1; % 0 for diagonal matrix, 1 for full matrix.
max_lags = 2; % number of lags to be considered for serial correlation of measurement errors
n_lag = 0; % Assume no AR in the first iteration.
detrend_price = "no"; % "yes" for deseasonalisation, otherwise "no"

par_names_model2 = define_parameters(LT, ncontracts, correlation, n_lag)';
n_par = size(par_names_model2, 1);
[par_optim_model2, log_L_optim_model2, par_init, trend, season, att, ytt, ett] = param_estim(y, ttm, deltat, detrend_price, n_par, ncontracts, par_names_model2, LT, correlation, max_lags);

n_lag = max_lags;
par_names_model3 = define_parameters(LT, ncontracts, correlation, n_lag)';
n_par = length(par_names_model3);
par_init = par_optim_model2;
[par_optim_model3, log_L_optim_model3, par_init, trend, season, att, ytt, ett] = param_estim_arp(y, ttm, deltat, detrend_price, n_par, par_names_model3, LT, correlation, par_init);

model2_result = set_parameters(LT, ncontracts, par_names_model2, par_model2, correlation, "yes");
model3_result = set_parameters(LT, ncontracts, par_names_model3, par_model3, correlation, "yes");

AIC_model2 = 2*log_L_model2 + par_optim(length(par_names_model2));
AIC_model3 = 2*log_L_model3 + par_optim(length(par_names_model3));

