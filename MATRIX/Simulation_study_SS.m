%% Simulation study of Schwartz-Smith Two-Factor model with serially correlated measurement errors

% This code runs the simulation study of the two-factor model, with AR(p)
% measurement errors.
% i.e. y_t = d_t + B_t * x_t + v_t,
% v_t = phi_1 * v_{t-1} + phi_2 * v_{t-2} + ... + phi_p * v_{t-p} + epsilon_t,
% where phi_p are diagonal matrices of AR coefficients at order p.
close all; clear all;

%% Simulation Study

% To simulate the data, use below codes. Otherwise, move to Line 65.

% Setting relevant inputs
T = [(1:5)]; % Maturities in years, say 1 - 5 months
ncontracts = size(T,2); % Number of contracts to be simulated
nobsn = 1000; % Number of observations to be simulated
deltat = 1/260; % Time difference between each case.

% Model selection
correlation = 1;
% correlation = 0 for independent volatilities in measurement errors,
% correlation = 1 for inter-correlated volatilities in measurement errors.

LT = "GBM"; % = "GBM" if the long-term factor follows GBM, "OU" if mean-reverting.

% Model parameters (Set true values)
kappa = 0.2462;
sigmachi = 0.0748;
lambdachi = 0.00622;
if LT == "GBM"
    gamma = [];
elseif LT == "OU"
    gamma = 0.0001;
end
mu = 0.4973;
sigmaxi = 0.5077;
lambdaxi = 0.4527;
rho = 0.2960;
vsig = repelem(0.01, ncontracts); % Volatilities of measurement errors.
correl = repelem(0.8, ncontracts); % intercorrelations between measurement errors.
n_lags = 2; % AR order (p)
m = repelem(0.4, n_lags * ncontracts); % AR coefficients of vt. Ensure stationarity when setting coefficients.

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
[y, x, ttm, v] = simulateY(par, par_names, T, ncontracts, nobsn, deltat, LT, correlation, serial);

%% Using provided simulated data, (Remove line 67 if not using provided simulated dataset)

% load('y_1000_N10.mat'); % Choose any file that begins with "y".
[nobsn, ncontracts] = size(y);
deltat = 1/260; % Time difference between each data point.

% Define the parameters for estimation
LT = "OU"; % GBM or OU.
n_season = 0;
correlation = 1; % 0 for diagonal matrix, 1 for full matrix.
max_lags = 2; % number of lags to be considered for serial correlation of measurement errors
detrend_price = "no"; % "yes" for deseasonalisation, otherwise "no"

global save_ett save_ytt save_att save_vt

% For two-step approach: we start by asssuming no AR in the measurement
% error (Model 1)
par_names_model0 = define_parameters(LT, ncontracts, correlation, 0)';
n_par = length(par_names_model0);
[par_model1, log_L_model0, par_init, trend, season, att1, ytt1, ett1] = param_estim(y, ttm, deltat, detrend_price, n_par, par_names_model0, LT, correlation, max_lags);
par_model0 = par_model1(1:8+ncontracts);
AIC_model0 = 2*log_L_model0+ 2*(length(par_names_model0));

% Model 1 with AR fitted on residuals
serial = "yes";
par_names_model1 = define_parameters(LT, ncontracts, correlation, 1)';
log_L_model1 = kf_v2_arp(par_model1, par_names_model1, y, deltat, ttm, LT, correlation, serial);
AIC_model1 = 2*log_L_model1 + 2*(length(par_names_model1));

% Model 2, unified estimation.
n_lag = max_lags;
par_names_model2 = define_parameters(LT, ncontracts, correlation, n_lag)';
n_par = length(par_names_model2);
par_init = par_model1;
[par_model2, log_L_model2, par_init, trend, season, att2, ytt2, ett2] = param_estim_arp(y, ttm, deltat, detrend_price, n_par, par_names_model2, LT, correlation, par_init);
AIC_model2 = 2*log_L_model2 + 2*(length(par_names_model2));

model0_result = set_parameters(LT, ncontracts, par_names_model0, par_model0, correlation, "no");
model1_result = set_parameters(LT, ncontracts, par_names_model1, par_model1, correlation, "yes");
model2_result = set_parameters(LT, ncontracts, par_names_model2, par_model2, correlation, "yes");