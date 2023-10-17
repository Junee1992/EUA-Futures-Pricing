%% Simulation study of Schwartz-Smith Two-Factor model with serially correlated measurement errors

% This code runs the simulation study of the two-factor model, with AR(p)
% measurement errors.
% i.e. y_t = d_t + B_t * x_t + v_t,
% v_t = phi_1 * v_{t-1} + phi_2 * v_{t-2} + ... + phi_p * v_{t-p} + epsilon_t,
% where phi_p are diagonal matrices of AR coefficients at order p.
close all; clear all;

%% Using provided simulated data, 

% load('data.mat'); % Choose any file that begins with "y".
[nobsn, ncontracts] = size(y);
deltat = 1/260; % Time difference between each data point.

% Define the parameters for estimation
LT = "OU"; % GBM or OU.
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