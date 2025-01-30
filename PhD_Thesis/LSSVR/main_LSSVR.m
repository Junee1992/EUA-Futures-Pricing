%% MLSSVR Model
% This code runs the Multi-output Least Sqaures Support Vector Regression
% Model. Assume futures contracts mature annually.

close all; clear all;
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

horiz = 20; 
output = forecast_LSSVR(y, horiz);

% Forecasting Plots
plot(nobsn-2*horiz+1:nobsn, y(nobsn-2*horiz+1:nobsn,1), 'k');
hold on
plot(nobsn-horiz:nobsn, output.y_data_temp(nobsn-horiz:nobsn,1), 'r')
plot(nobsn-horiz:nobsn, [y(nobsn-horiz,1); output.y_data_temp(nobsn-horiz+1:nobsn,1) + 1.96 * sqrt(output.res_lssvr_var(:,1))], 'r--')
plot(nobsn-horiz:nobsn, [y(nobsn-horiz,1); output.y_data_temp(nobsn-horiz+1:nobsn,1) - 1.96 * sqrt(output.res_lssvr_var(:,1))], 'r--')