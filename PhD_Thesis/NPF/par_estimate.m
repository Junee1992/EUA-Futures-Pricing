%% Parameter Estimation
function [par_optim, log_L_optim, par_init, trend, season, att, ytt, ett, ett_1, Ptt, fitresult_linear, fitresult_season] = par_estimate(y, ttm, model_options, max_lags)

fields = fieldnames(model_options);
for i = 1:numel(fields)
    fieldName = fields{i};
    eval([fieldName ' = model_options.(fieldName);']);
end

[nobsn, ncontracts] = size(y);
% Initial estimation - assume no AR process
n_lag = 0; n_par = 7 + (LT == "OU") + ncontracts + ncontracts * (correlation == 1);
serial = "no";

% Parameter constraints
[A, b, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, n_lag);
% lb(4) = 0.0001; lb(1) = 0.0001;
if detrend_price == "yes"
    % Detrend the price
    [y_detrend, trend, fitresult_linear] = linear_detrend(y);
    % [season, y_deseason, fitresult_season, gof] = deseasonalise_v2(y_detrend);
    seasonality = "sinusoid_b";
    for i = 1:10
        [season, p, fitresult_season, gof] = deseasonalise(y_detrend, i, seasonality);
        rmse_seasonality(i) = gof.rmse;
    end
    no_season = min(find(rmse_seasonality == min(rmse_seasonality)));
    [season, p, fitresult_season, gof] = deseasonalise(y_detrend, no_season, seasonality);
    y_deseason = y_detrend - season;
else
    trend = []; season = []; fitresult_linear = []; fitresult_season = [];
    y_deseason = y;
end

% Grid search
correl = correlation;
init = grid_search(lb, ub, n_lag, ncontracts, LT, correl);
init_idx = find(init(:,1) < init(:,4));
init(init_idx,:) = [];
% log_L = kf_v1(init(1,:), par_names, y_deseason, deltat, ttm, LT, correl, serial);
% parpool(8);
parfor i = 1:size(init,1)
    log_L = kf_SS(init(i,:), y_deseason, ttm, model_options);
    logL(i) = log_L;
end
par_init = mean(init(find(logL == min(logL)),:),1);

for i = 1:length(par_init)
    if lb(i) >= par_init(i)
        par_init(i) = par_init(i) + 1e-4;
    elseif par_init(i) >= ub(i)
        par_init(i) = par_init(i) - 1e-4;
    end
end

if correlation == 1
    rho_range = linspace(-0.9999, 0.9999, 5);
    [grids{1:ncontracts}] = ndgrid(rho_range);
    parameter_grid = cell2mat(cellfun(@(x) x(:), grids, 'UniformOutput', false));
    init = [repmat(par_init(:,1:8+ncontracts),length(parameter_grid), 1) parameter_grid];
    parfor i = 1:length(init)
        log_L = kf_SS(init(i,:), y_deseason, ttm, model_options);
        logL(i) = log_L;
    end
    % par_init = mean(init(find(logL == min(logL)),:),1);
    par_init = init(find(logL == min(logL)),:);
    par_init = par_init(1,:);
end

for i = 1:length(par_init)
    if lb(i) >= par_init(i)
        par_init(i) = par_init(i) + 1e-4;
    elseif par_init(i) >= ub(i)
        par_init(i) = par_init(i) - 1e-4;
    end
end

global save_ett save_ytt save_att save_vt save_att_1 save_Ptt save_Ptt_1 save_xt_1n save_Pt_1n save_ett_1

% Setting options for optimisation
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 500, 'MaxFunEvals', 50000);

if LT == "GBM"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_SS, par_init, [], [], [], [], lb, ub, [], options, y_deseason, ttm, model_options);
elseif LT == "OU"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_SS, par_init, A, b, [], [], lb, ub, [], options, y_deseason, ttm, model_options);
end

% Checking AR
if max_lags > 0
    for j = 1:ncontracts
        ARmodel = estimate(arima('Constant', 0, 'ARLags', [], 'Distribution', 'Gaussian'), save_ett(:,j), 'display', 'off');
        ARmodelAIC(1, j) = ARmodel.summarize.AIC;
        for l = 1:max_lags
            try
            [ARmodel] = estimate(arima('Constant', 0, 'ARLags', 1:l, 'Distribution', 'Gaussian'), save_ett(:,j), 'display', 'off');
            catch
                  [ARmodel] = estimate(arima('Constant', 0, 'ARLags', 1:1, 'Distribution', 'Gaussian'), save_ett(:,j), 'display', 'off');
            end
            AR_coef{l+1,j} = ARmodel.AR;
            ARmodelAIC(l+1,j) = ARmodel.summarize.AIC;
        end
    end
    contractNames = arrayfun(@(x) sprintf('C%d', x), 1:ncontracts, 'UniformOutput', false);
    AICsummary = array2table([[0:max_lags]', ARmodelAIC], 'VariableNames', ['Lags' contractNames]);
    for j = 1:ncontracts
        lags_temp(j) = min(find(ARmodelAIC(:,j) == min(ARmodelAIC(:,j))))-1; 
    end

    n_lag = max(lags_temp); % number of lags for serial correlation of measurement errors
    m_init = zeros(ncontracts, n_lag);
    for j = 1:ncontracts
        m_init(j,:) = [cell2mat(AR_coef{lags_temp(j)+1,j}), repelem(0, n_lag-lags_temp(j))];
    end
    par_optimised = [par_optimised reshape(m_init, 1, size(m_init,1)*size(m_init,2))];
else
    par_optimised = par_optimised;
end

par_optim = par_optimised;
log_L_optim = log_L;
att = save_att;
ytt = save_ytt;
ett = save_ett;
vtt = save_vt;
ett_1 = save_ett_1;
Ptt = save_Ptt;

% if LT == "GBM"
%     G = [exp(-par_optim(1)* deltat), 0; 0, 1];
% elseif LT == "OU"
%     G = [exp(-par_optim(1) * deltat), 0; 0, exp(-par_optim(4)* deltat)];
% end
% 
% if LT == "GBM"
%     a0 = [0; mean(y_deseason(1,:), 2)];
%     P0 = [0.01, 0.01; 0.01, 0.01];
% 
% elseif LT == "OU"
%     a0 = [0; mean(y_deseason(1,:), 2)];
%     P0 = [0.01, 0.01; 0.01, 0.01];
% 
% end
% [save_xt_1n, save_Pt_1n] = kalman_sm(save_Ptt, save_Ptt_1, save_att, save_att_1, G, a0, P0);

end