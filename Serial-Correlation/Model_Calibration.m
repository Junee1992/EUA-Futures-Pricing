%% Calibration of EUA Futures Prices

% This code runs the model calibration of deseasonalised EUA Futures
% Prices.

clear all; close all;

% Loading the deseasonalised prices
load Deseasonalised_price_EUA_Dec2021.mat
y = detrend_price;
[nobsn, ncontracts] = size(y);
deltat = 1/360;
% nobsn - number of observations.
% ncontracts - number of contracts.
% ttm - time to maturity.
% deltat - the time difference between each observation, assumed to be
% 1/360 for daily price data.

% Select your model
model = "price";
% Select model = "price" for modelling of logarithmic price.
% Select model = "return" for modelling of logarithmic return.

% Model selection
ms_error = 2;
% ms_error = 1 for no independent volatilities in measurement errors,
% ms_error = 2 for inter-correlated volatilities in measurement errors.

% Choice of having an AR process
AR_process = "noAR";
% Selected AR_process = "AR" to consider AR process in measurement errors,
% and for estimation of AR coefficients.
% Otherwise, use AR_process = "noAR" for serially uncorrelated measurement errors.

% Initial values for parameters
kappa_init = 1.5;
sigmachi_init = 0.05;
lambdachi_init = 0.02;
gamma_init = 0.5;
mu_init = 1;
sigmaxi_init = 0.05;
lambdaxi_init = 0.02;
rho_init = 0.9;
vsig_init = repelem(0.1, ncontracts); % Volatilities of measurement errors.
m_init = repelem(-0.1, ncontracts); % AR coefficients
correl_init = repelem(0.9, ncontracts); % Intercorrelations

% Parameter estimation contraint to avoid parameter identification problem
% Modify boundaries accordingly for each model
A = [-1, 0, 0, 1, 0, 0, 0, 0, repelem(0, ncontracts * 2)];
b = 0;
par_init = [kappa_init, sigmachi_init, lambdachi_init, gamma_init, mu_init, sigmaxi_init, lambdaxi_init, rho_init, vsig_init, correl_init];

% Setting lower and upper bounds for each parameter
lb(1:8) = [0.0001, 0.0001, 0.0001, 0.0001, -4, 0, 0, -0.9999];
lb(9:8+ncontracts) = 0.000001;
lb(9+ncontracts:8+2*ncontracts) = -0.9999;
ub(1:8) = [5, 5, 0.9999, 5, 5, 5, 0.9999, 0.9999];
ub(9:8+ncontracts) = 1;
ub(9+ncontracts:8+2*ncontracts) = 0.9999;

% Initial distribution of state variables

if model == "price"
    a0 = [0 ; mu_init / gamma_init]; % a_1|0, initial expectation
    P0 = [sigmachi_init^2 / (2*kappa_init), sigmachi_init * sigmaxi_init * rho_init / (kappa_init + gamma_init);
        sigmachi_init * sigmaxi_init * rho_init / (kappa_init + gamma_init), sigmaxi_init^2 / (2 * gamma_init)];% P_1_0, initial variance matrix
    
elseif model == "return"
    a0 = [0; y(2,1); 0; y(1,1)];
    % a_1|0, initial expectation, assumed to take the first two prices from
    % first futures contract.
    P0 = [0.1, 0.1, 0, 0; 0.1, 0.1, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0];
    % P_1|0, initial variance matrix, assumed to be in the matrix
    % form similar to w^r in equation (14) of the article.
end

% Setting options for optimisation
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 1000, 'MaxFunEvals', 50000);

% Optimisation
global save_att save_ytt save_ett
[par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@KFS_ar, par_init, A, b, [], [], lb, ub, [], options, y, deltat, ttm, model, a0, P0, ms_error, AR_process);

% Checking for AR(1) process in the marginal measurement error
[m, AIC0, AIC1] = AR1_fit(save_ett);

if norm(m) > 0
    AR_process == "AR"
    par_optimised = [par_optimised, m];
    
    % Re-setting the lower and upper bounds of parameter estimates
    for i = 1:ncontracts
        if m_init(i) == 0
            lb_m = 0;
            ub_m = 0;
            lb = [lb lb_m];
            ub = [ub ub_m];
        else
            lb_m = -0.9999;
            ub_m = 0.9999;
            lb = [lb lb_m];
            ub = [ub ub_m];
        end
    end
    
    A = [A, repelem(0, ncontracts)];
    [par_optimised_ar, log_L_ar, exitflag, output, lambda, grad, hessian] = fmincon(@KFS_ar, par_optimised, A, b, [], [], lb, ub, [], options, y, deltat, ttm, model, a0, P0, ms_error, AR_process);
else
    par_optimised = par_optimised;
end

% Results

par_names = ["kappa", "sigmachi", "lambdachi", "gamma", "muxi", "sigmaxi", "lambdaxi", "rho"]';
par_table = table(par_names, par_optimised(1:8)')

vol_names = string(1:ncontracts)';
vol_table = table(vol_names, par_optimised(9:8+ncontracts)')

correl_names = string(1:ncontracts)';
correl_table = table(correl_names, par_optimised(9+ncontracts:8+2*ncontracts)')

if AR_process == "AR"
    AR_names = string(1:ncontracts)';
    AR_table = table(AR_names, par_optimised(9+2*ncontracts:8+3*ncontracts)')
end
