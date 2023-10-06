## This file includes all instructions on how to use codes provided in this Github repository.

## define_parameters.m

par_names = define_parameters(LT, ncontracts, correlations, n_lags)

Lists the name of parameters to be estimated, according to the setting selected by the user.
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- ncontracts -> the number of available contracts in the price dataset.
- correlation -> 0, assumes measurement errors are independent
              -> 1, assumes measurement errors are intercorrelated
- n_lag -> the number of lags to be considered for estimation of serial correlations of errors.

## deseasonalise.m

[seasonal, periodicity, fitresult, gof] = deseasonalise(price, n_season, seasonality)

Performs deseasonalisation of the data via examining the power spectral density of the mean of data.

Input:
- price -> the price data
- n_season -> the number of seasonal components to be considered
- seasonality -> "sinusoid_a", using f(x) = SUM {a_t * cos(t*2*pi*x) + b_t * sin(t*2*pi*x)} for t = 1, ..., n_season.
              -> "sinusoid_b", using f(x) = SUM {a_t * cos({periodicity}*2*pi*x) + b_t * sin({periodicity}*2*pi*x)}, where the periodicity is determined
                 by the spectral density

Output:
- season -> seasonal component
- periodicity -> periodicity determined by the spectral density
- fitresult -> the model of the seasonal function
- gof -> results of the goodness of fit.

## grid_search.m

init = grid_search(lb, ub, n_lag, ncontracts, LT, correlation)

Performs the grid search procedure on key model parameters (kappa, sigmachi, gamma, mu, sigmaxi, rho). This procedure assumes three grid points
for each parameter.

Input:
- lb -> lower bound of parameter space
- ub -> upper bound of parameter space
- n_lag -> the number of lags to be considered for estimation of serial correlations of errors.
- ncontracts -> the number of available contracts in the price dataset
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- correlation -> 0, assumes measurement errors are independent
              -> 1, assumes measurement errors are intercorrelated

## kf_v1.m

log_L = kf_v1(par_init, par_names, y, deltat, ttm, LT, correlation, serial)

Performs Kalman filtering and computes the log-likelihood function, assuming serial independence in measurement errors.

Input:
- par_init -> initial set of parameters
- par_names -> the list of names of model parameters
- y -> the price data
- deltat -> the difference between each price. Set deltat = 1/260 for daily data.
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- correlation -> 0, assumes measurement errors are independent
              -> 1, assumes measurement errors are intercorrelated
- serial -> "no", assumes serial independence in measurement errors
         -> "yes", assumes serial correlation in measurement errors

Output:
- log_L -> the value of the log-likelihood function

## kf_v3_arp.m

log_L = kf_v3_arp(par_init, par_names, y, deltat, ttm, LT, correlation, serial)

Performs Kalman filtering and computes the log-likelihood function, assuming serial dependence in measurement errors.

Input:
- par_init -> initial set of parameters
- par_names -> the list of names of model parameters
- y -> the price data
- deltat -> the time difference between each price. Set deltat = 1/260 for daily data.
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- correlation -> 0, assumes measurement errors are independent
              -> 1, assumes measurement errors are intercorrelated
- serial -> "no", assumes serial independence in measurement errors
         -> "yes", assumes serial correlation in measurement errors

Output:
- log_L -> the value of the log-likelihood function

## ks_v2.m

wt_1 = ks_v2(n_t, k, C, G, att, att_1, Ptt, Ptt_1)

Performs Kalman smoothing that is required for construction of likelihood function in case when p > 2.

Input: 
- n_t -> the current final time
- k -> lags considered in a particular iteration
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- par -> a set of parameters
- deltat -> the time difference between each price. Set deltat = 1/260 for daily data.
- ttm -> time to maturity for all prices at each time point
- att -> E[X_t|F_t], updated state
- att_1 -> E[X_t|F_{t-1}], state prediction
- Ptt -> Var[X_t|F_t], updated state variance
- Ptt_1 -> Var[X_t|F_{t-1}], state variance prediction

Output:
- wt_1 -> E[w_{t-k+1}|F_{t-1}], an estimate of the state error

## linear_detrend.m

[y_detrend, trend, fitresult] = linear_detrend(y)

Performs de-linearisation of the data, based on the linear regression.

Input:
- y -> the data.

Output:
- y_detrend -> the detrended data
- trend -> the trend
- fitresult -> the model of the linear trend

## param_constr.m

[A, b, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, n_lag)

Sets constraints for model parameters.

Input:
- n_par -> the number of parameters
- ncontracts -> the number of available contracts in the price dataset
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- correlation -> 0, assumes measurement errors are independent
              -> 1, assumes measurement errors are intercorrelated
- n_lag -> the number of lags to be considered for estimation of serial correlations of errors

Output:
- A, b -> linear constraints, such that A*x <= b.
- lb, ub -> lower and upper bounds for each parameter

## param_estim.m

[par_optim, log_L_optim, par_init, trend, season, att, ytt, ett, fitresult_linear, fitresult_season] = param_estim(y, ttm, deltat, detrend_price, n_par, par_names, LT, correlation, max_lags)

Performs parameter estimation when measurement errors are serially independent.

Input:
- y -> the price data.
- ttm -> time to maturity for all prices at each time point
- deltat -> the time difference between each price. Set deltat = 1/260 for daily data.
- detrend_price -> "yes", if the user wants the data to be detrended.
                -> "no", if the user wants the original data to be considered.
- n_par -> the number of parameters
- par_names -> the list of names of model parameters
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- correlation -> 0, assumes measurement errors are independent
              -> 1, assumes measurement errors are intercorrelated
- max_lags -> the number of maximum lags to be examined on fitted measurement errors

Output:
- par_optim -> estimated parameters
- log_L_optim -> the minimum negative log likelihood function
- par_init -> initial parameters used after the grid search
- trend -> the trend
- season -> seasonal component
- att -> the estimated state variables
- ytt -> the estimated price
- ett -> the fitted residuals
- fitresult_linear -> the model of the linear function
- fitresult_season -> the model of the seasonal function

## param_estim_arp.m

[par_optim, log_L_optim, par_init, trend, season, att, ytt, ett, vtt, fitresult_linear, fitresult_season] = param_estim_arp(y, ttm, deltat, detrend_price, n_par, par_names, LT, correlation, par_init)

Performs parameter estimation when measurement errors are serially correlated.

Input:
- y -> the price data.
- ttm -> time to maturity for all prices at each time point
- deltat -> the time difference between each price. Set deltat = 1/260 for daily data.
- detrend_price -> "yes", if the user wants the data to be detrended.
                -> "no", if the user wants the original data to be considered.
- n_par -> the number of parameters
- par_names -> the list of names of model parameters
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- correlation -> 0, assumes measurement errors are independent
              -> 1, assumes measurement errors are intercorrelated
- par_init -> initial set of parameters

Output:
- par_optim -> estimated parameters
- log_L_optim -> the minimum negative log likelihood function
- par_init -> initial parameters used after the grid search
- trend -> the trend
- season -> seasonal component
- att -> the estimated state variables
- ytt -> the estimated price
- ett -> the fitted residuals
- vtt -> E[v_t|F_t], the measurement errors that follows AR process
- fitresult_linear -> the model of the linear function
- fitresult_season -> the model of the seasonal function

## set_parameters.m

param = set_parameters(LT, ncontracts, par_names, par_init, correlation, serial)

Tabularises the values of parameters with the names of parameters.

Input:
- LT -> "GBM" if the long-term factor follows geometric Brownian motion.
     -> "OU" if the long-term factor follows Ornstein-Uhlenbeck process.
- ncontracts -> the number of available contracts in the price dataset
- par_names -> the list of names of model parameters
- par -> set of parameters
- correlation -> 0, assumes measurement errors are independent
              -> 1, assumes measurement errors are intercorrelated
- serial -> "no", assumes serial independence in measurement errors
         -> "yes", assumes serial correlation in measurement errors
