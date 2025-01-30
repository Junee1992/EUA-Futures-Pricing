%% NPF (Nested Particle Filter) Model
% This code runs the NPF (Nested Particle Filter) algorithm for estimation
% of states and parameters. Assume futures contracts mature annually.

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
err = "laplace";
% = "laplace" for Multivariate Generalised Laplace errors,
% = "hyperbolic" for Multivariate Generalised Hyperbolic errors
if err == "laplace"
    alpha = 2;
    beta = alpha;
elseif err == "hyperbolic"
    lambda = 0.5;
    psi = 2;
    chi = 0.8;
    nu = repelem(0.005, ncontracts);
    beta = [lambda chi psi nu];
else
    beta = [];
end

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
if err == "laplace"
    par_names = [par_names "sG"];
elseif err == "hyperbolic"
    nu_name = [];
    for i = 1:ncontracts
        nu_name = [nu_name sprintf("nu_%d", i)];
    end
    par_names = [par_names "lambda" "psi" "chi" nu_name];
end

par = [kappa sigmachi lambdachi gamma mu sigmaxi lambdaxi rho_chixi s rho phi beta];
% par = set_parameters(LT, ncontracts, par_names, param, correlation, serial, err);
T = 1:ncontracts; detrend_price = "no";
model_options = struct('LT', LT, 'correlation', correlation, 'nobsn', nobsn, ...
    'ncontracts', ncontracts, 'par_names', par_names, 'deltat', deltat, 'T', T, ...
    'detrend_price', detrend_price, 'err', err);

% Simulate the data.
[y, x, ttm, v] = simulateY_nG(par, model_options);

% Forecast
n_forecast = 20; M = 1000; N = 1000;
max_lags = 0;
model_options = struct('LT', LT, 'correlation', correlation, 'nobsn', nobsn, ...
    'ncontracts', ncontracts, 'par_names', par_names, 'deltat', deltat, 'T', T, ...
    'detrend_price', detrend_price, 'err', err, 'n_forecast', n_forecast, ...
    'M', M, 'N', N, 'max_lags', max_lags);
output = forecast_NPF(y, ttm, model_options);

% Plot of forecast (1st available contract)
figure;
y_temp = output.y_temp; varn = output.varn;
plot(y(:,1), 'k');
hold on
plot(nobsn-n_forecast:nobsn, y_temp(nobsn-n_forecast:nobsn,1), 'r');
for i = 1:n_forecast
    lb_y(i,:) = y_temp(nobsn-n_forecast+i,:)' - 1.96 * sqrt(diag(varn(:,:,i)));
    ub_y(i,:) = y_temp(nobsn-n_forecast+i,:)' + 1.96 * sqrt(diag(varn(:,:,i)));
end
plot(nobsn-n_forecast:nobsn, [y_temp(nobsn-n_forecast,1); lb_y(:,1)], 'r--');
plot(nobsn-n_forecast:nobsn, [y_temp(nobsn-n_forecast,1); ub_y(:,1)], 'r--');

