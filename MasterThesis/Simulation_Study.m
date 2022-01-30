%% Code for Simulation Study
% On Multi-Factor State-Space Modelling and Forecasting of EUA Futures
% Prices (Han, 2020)

% This file contains the code used for Jun Han's {Master of Research}
% thesis, submitted to Macquarie University, Australia, in Nov 2020.

% The code provides parameter estimation of a simulated datset based on
% modified Schwartz-Smith two-factor model, allowing for correlated
% measurement errors. Other options are also provided.

% This code is divided into four different sections.

% First section requires the user to put required inputs, such as parameters,
% number of observations, number of contracts, maturities, and the type of
% the model for consideration.

% Second section performs formulating possible initial values for grid
% search. The user is required to type the number of grid points.

% Third section performs parameter estimation, using possible initial
% values from second section. The user is recommended for use of parallel
% computing in this section.

% Fourth section computes the standard error of parameters based on the
% best initial and estimated values of parameters. The user need to request
% the number of simulations for initialising the bootstrapping approach in
% computing those standard errors. The user is recommended for use of
% parallel computing in this section.

%% Required Inputs
clear all; clc; close all;

model_selection = 3; % Choose model_selection = 1 for GBM model,
%        model_selection = 2 for SSM model,
%        model_selection = 3 for extended SSM model.

vary_ME = 1; % Choose vary_ME = 1 for no correlations between measurement errors,
%        vary_ME = 2 for estimating volatilities as a
%        parametric function,
%        vary_ME = 3 for allowing correlations between
%        measurement errors,
%        vary_ME = 4 for allowing correlations between
%        measurement errors, and estimating volatilities as a
%        parametric function.

yeardays = 360; % Number of days in one year.
monthdays = 30; % Number of days in a month.
deltat = 1/yeardays; % Time difference for each observation.
T = [1:12]/12; % Maturities in years.
nobsn = 10; % Number of observations per contract.
ncontracts = size(T,2); % Number of contracts, must be consistent with number of maturities considered.

% Required parameters according to state and measurement equations
kappa = 2;
sigmachi = 0.2;
lambdachi = 0.02;
gamma = 1;
mu = 0.5;
sigmaxi = 0.3;
lambdaxi = 0.03;
rho = -0.5;
param = [kappa, sigmachi, lambdachi, gamma, mu, sigmaxi, lambdaxi, rho];

% Volatilities and correlations of measurement errors

% Volatilities
s = [repelem(0.01, ncontracts)];

% Parametrised volatilities
% s = v1 + v2 * exp(v3 * T), where T is a vector of maturities in years.
v1 = 0.0001;
v2 = 0.0001;
v3 = 5;
v = [v1, v2, v3];

% Correlations
correl = [repelem(0.9, ncontracts)];

% Setting the parameter set
if vary_ME == 1
    par = [param, s];
    
elseif vary_ME == 2
    par = [param, v];
    
elseif vary_ME == 3
    par = [param, s, correl];
    
elseif vary_ME == 4
    par = [param, v, correl];
    
else
    error("Please choose the correct vary_ME option")
end

% If time-to-maturity is calculated by the user, use
% ttm_by_user = ttm; % ttm is user's time-to-maturity according to the
% known trading dates and expiry dates.
ttm_by_user = 0; % If the user wants to only simulate the dataset without any given information.

% Simulate a dataset
[simY, x, ttm] = simulatePrice(par, T, yeardays, monthdays, ncontracts, nobsn, deltat, model_selection, vary_ME, ttm_by_user);
% simY: the simulated futures price based on the state and the
% measurement equations.
% x: the values of latent variables
% ttm: the time-to-maturity according to the maturity times given in the
% input.

% Setting lower and upper bounds for parameters

% Main parameters
lb(1:8) = [0.0001, 0.0001, -5, 0.0001, -5, 0.0001, -5, -1];
ub(1:8) = [5, 5, 5, 5, 5, 5, 5, 1];

% Volatilities and correlations
if vary_ME == 1
    lb(9:8+ncontracts) = 0.000001;
    ub(9:8+ncontracts) = 1;
    
    if model_selection == 1
        lb([1 2 3 4 8]) = 0; ub([1 2 3 4 8]) = 0;
        
    elseif model_selection == 2
        lb(4) = 0;
        ub(4) = 0;
        
    elseif model_selection == 3
        %Parameter Estimation Constraint
        A = [-1, 0, 0, 1, 0, 0, 0, 0, repelem(0, ncontracts)];
        b = 0;
        
    else
        print('Please choose the correct option')
    end
    
elseif vary_ME == 2
    lb(9:10) = 0.000001;
    lb(11) = -10;
    ub(9:10) = 1;
    ub(11) = 10;
    
    if model_selection == 1
        lb([1 2 3 4 8]) = 0; ub([1 2 3 4 8]) = 0;
        
    elseif model_selection == 2
        lb(4) = 0; ub(4) = 0;
        
    elseif model_selection == 3
        A = [-1, 0, 0, 1, 0, 0, 0, 0, repelem(0, 3)];
        b = 0;
        
    else
        print('Please choose the correct option')
    end
    
elseif vary_ME == 3
    lb(9:8+ncontracts) = 0.000001;
    lb(9+ncontracts:8+2*ncontracts) = -1;
    ub(9:8+ncontracts) = 1;
    ub(9+ncontracts:8+2*ncontracts) = 1;
    
    if model_selection == 1
        lb([1 2 3 4 8]) = 0; ub([1 2 3 4 8]) = 0;
        
    elseif model_selection == 2
        lb(4) = 0; ub(4) = 0;
        
    elseif model_selection == 3
        A = [-1, 0, 0, 1, 0, 0, 0, 0, repelem(0, 2*ncontracts)];
        b = 0;
        
    else
        print('Please choose the correct option')
    end
    
elseif vary_ME == 4
    lb(9:10) = 0.000001;
    lb(11) = -10;
    lb(12:11+ncontracts) = -1;
    ub(9:10) = 1;
    ub(11) = 10;
    ub(12:11+ncontracts) = 1;
    
    if model_selection == 1
        lb([1 2 3 4 8]) = 0; ub([1 2 3 4 8]) = 0;
        
    elseif model_selection == 2
        lb(4) = 0; ub(4) = 0;
        
    elseif model_selection == 3
        A = [-1, 0, 0, 1, 0, 0, 0, 0, repelem(0, 3+ncontracts)];
        b = 0;
        
    else
        print('Please choose the correct option')
    end
end


%% Grid search
ngrid = 2; % Choose the number of grid points.

mid = (lb + ub) / 2;
mp1 = (mid + lb) / 2;
mp2 = (mid + ub) / 2;

if model_selection == 1
    npar = 3;
elseif model_selection == 2
    npar = 7;
elseif model_selection == 3
    npar = 8;
end

if ngrid == 2
    grid = [mp1; mp2]';
elseif ngrid == 3
    grid = [mp1; mid; mp2]';
end

est = zeros(ngrid^npar, length(par)+1);
init = [];

if model_selection == 1
    for m = grid(5,:)
        for sx = grid(6,:)
            for lx = grid(7,:)
                init = [init; 0, 0, 0, 0, m, sx, lx, 0, repelem(0.01, ncontracts)];
            end
        end
    end
    
elseif model_selection == 2
    for k = grid(1,:)
        for sc = grid(2,:)
            for lc = grid(3,:)
                for m = grid(5,:)
                    for sx = grid(6,:)
                        for lx = grid(7,:)
                            for rh = grid(8,:)
                                init = [init; k, sc, lc, 0, m, sx, lx, rh, repelem(0.01, ncontracts)];
                            end
                        end
                    end
                end
            end
        end
    end
    
elseif model_selection == 3
    for k = grid(1,:)
        for sc = grid(2,:)
            for lc = grid(3,:)
                for g = grid(4,:)
                    for m = grid(5,:)
                        for sx = grid(6,:)
                            for lx = grid(7,:)
                                for rh = grid(8,:)
                                    init = [init; k, sc, lc, g, m, sx, lx, rh, repelem(0.01, ncontracts)];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if vary_ME == 2
    init = [init(:,1:8), repmat(0.001, size(init,1), 2), repmat(5, size(init,1),1)];
    
elseif vary_ME == 3
    init = [init, repmat(0.9, size(init,1), ncontracts)];
    
elseif vary_ME == 4
    init = [init(:,1:8), repmat(0.001, size(init,1), 2), repmat(5, size(init,1),1), repmat(0.9, size(init,1), ncontracts)];

end

%% Parameter estimation
options = optimset('Display', 'iter','TolFun',1e-06,'TolX',1e-06 ,'MaxIter',1000,'MaxFunEvals',4000);

for i = 1:ngrid^npar
    i
    par0 = init(i,:);
    if any(model_selection == [1,2])
        [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@KFS, par0, [], [], [], [], lb, ub, [], options, simY, deltat, ttm, model_selection, vary_ME, T);
        est(i,:) = [par_optimised, log_L];
    elseif model_selection == 3
        [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@KFS, par0, A, b, [], [], lb, ub, [], options, simY, deltat, ttm, model_selection, vary_ME, T);
        est(i,:) = [par_optimised, log_L];    
    end
    % Note, for GBM model, RTS smoother does not work due to prediction
    % covariance being singular.
end
disp('RTS smoother not produced for Model 1')

index = (est(:,end) == min(est(:,end)));
par_est = est(index,:); % Parameter estimates based on the best initial values, along with the log likelihood.
best_init = init(index,:); % Best initial values

global save_xf save_xs_MBF save_xs_RTS save_ett save_ytt

% save_xf - estimated values of latent variables using Kalman filter
% save_xs_RTS - estimated values of latent variables using RTS smoother
% save_xs_MBF - estimated values of latent variables using MBF smoother
% save_ett - residuals, calculated as (actual price - fitted price)
% save_ytt - fitted values of futures prices
% To obtain these estimates, the user needs to re-run fmincon using the
% best initial value.

if any(model_selection == [1,2])
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@KFS, best_init, [], [], [], [], lb, ub, [], options, simY, deltat, ttm, model_selection, vary_ME, T);
elseif model_selection == 3
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@KFS, best_init, A, b, [], [], lb, ub, [], options, simY, deltat, ttm, model_selection, vary_ME, T);
end

%% Standard errors of parameter estimates
h = 2000; % Number of resampling for computing standard errors
save_par = zeros(h,size(best_init,2));
    
for i = 1:h
    i
    simY1 = simulatePrice(best_est(1:end-1), T, yeardays, monthdays, ncontracts, nobsn, deltat, model_selection, vary_ME, ttm);
    if any(model_selection == [1,2])
        [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@KFS, best_init, [], [], [], [], lb, ub, [], options, simY1, deltat, ttm, model_selection, vary_ME, T);
    elseif model_selection == 3
        [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@KFS, best_init, A, b, [], [], lb, ub, [], options, simY1, deltat, ttm, model_selection, vary_ME, T);
    end
    save_par1(i,:) = par_optimised;
end

std_error = sqrt(var(save_par1)/h);