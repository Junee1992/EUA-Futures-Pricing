%% SS Model - Kalman filter (with AR(p) measurement errors)
% This code runs the Schwartz-Smith Two-factor Model (2000) with two
% factors following correlated bivariate Ornstein-Uhlenbeck processes.
% Assume futures contracts mature annually.

clear all; close all;

rng(1000);

% Model Options (Simulation)
LT = "OU"; correlation = 1; % = 1 if inter-correlated measurement errors, = 0 if diagonal.
deltat = 1/260; % Assume 260 days in a year.
ncontracts = 5; nobsn = 1000; % Number of contracts, number of observations, n x K
kappa = 0.5; sigmachi = 0.3; lambdachi = 0.02;
gamma = 0.1; mu = 0.5; sigmaxi = 0.2; lambdaxi = 0.02;
rho_chixi = 0.8;
s = repelem(0.02, ncontracts);
rho = repelem(0.9, ncontracts);
phi = [];

if ~isempty(rho)
    correlation = 1;
else
    correlation = 0;
end

if ~isempty(phi)
    n_lag = length(phi)/ncontracts;
    serial = "yes";
else
    n_lag = 0;
    serial = "no";
end

par_names = define_parameters(LT, ncontracts, correlation, n_lag);
param = [kappa sigmachi lambdachi gamma mu sigmaxi lambdaxi rho_chixi s rho phi];
par = set_parameters(LT, ncontracts, par_names, param, correlation, serial, "normal");
T = 1:ncontracts; detrend_price = "no";
model_options = struct('LT', LT, 'correlation', correlation, 'nobsn', nobsn, ...
    'ncontracts', ncontracts, 'par_names', par_names, 'deltat', deltat, 'T', T, ...
    'detrend_price', detrend_price);

% Simulate the data.
[y, x, ttm, v] = simulate_Y(par, model_options);

% Global variables
global save_ett save_ytt save_att save_vtt

% Data Calibration / Model Estimation
% att: estimated states
% ytt: estimated measurements
% ett: measurement errors
% par_optim: optimised parameters
% max_lags = 0; % Maximum lags to be considered initially.
% [par_optim, ~, ~, ~, ~, att, ytt, ett] = par_estimate(y, ttm, model_options, max_lags);

% Forecasting
n_forecast = 20;
y_pred = [];
n_temp = nobsn - n_forecast;
y_temp = y(1:n_temp,:); y_fore = [];
detrend_price = "yes";
LT = "OU"; % GBM or OU.
n_season = 0;
correlation = 1; % 0 for diagonal matrix, 1 for full matrix.
max_lags = 1; % number of lags to be considered for serial correlation of measurement errors
n_lag = 0; % Assume no AR in the first iteration.
par_names = define_parameters(LT, ncontracts, correlation, n_lag)';
model_options = struct('LT', LT, 'correlation', correlation, 'par_names', par_names, 'deltat', deltat, ...
    'detrend_price', detrend_price, 'max_lags', max_lags, 'n_forecast', n_forecast, 'n_temp', n_temp);
output = forecast_KF(y_temp, ttm, model_options);
