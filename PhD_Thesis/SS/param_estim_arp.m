%% Parameter Estimation
function [par_optim, log_L_optim, par_init, trend, season, att, ytt, ett, vtt, fitresult_linear, fitresult_season, att_1, Ptt, Ptt_1] = param_estim_arp(y, ttm, deltat, detrend_price, n_par, par_names, LT, correlation, par_init)

[nobsn, ncontracts] = size(y);
if sum(contains(par_names, "phi")) > 0
    n_lag = sum(contains(par_names, "phi"))/ncontracts;
    serial = "yes";
else
    n_lag = 0; serial = "no";
end

if LT == "GBM"
    model_par = 7;
elseif LT == "OU"
    model_par = 8;
end

% Parameter constraints
[A, b, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, n_lag);
idx = find(contains(par_names, "phi"));
for p = 1:length(idx)
    i = idx(p);
    if par_init(i) == 0
        lb(i) = 0; ub(i) = 0;
    end
    if par_init(i) > ub(i)
        par_init(i) = ub(i)-1e-4;
    end
end
% lb(9:8+ncontracts) = 1e-5;
if detrend_price == "yes"
    % Detrend the price
    [y_detrend, trend, fitresult_linear] = linear_detrend(y);
    seasonality = "sinusoid_b";
    for i = 1:10
        [season, p, fitresult_season, gof] = deseasonalise(y_detrend, i, seasonality);
        rmse_seasonality(i) = gof.rmse;
    end
    no_season = min(find(rmse_seasonality == min(rmse_seasonality)));
    [season, p, fitresult_season, gof] = deseasonalise(y_detrend, no_season, seasonality);
    y_deseason = y_detrend - season;
else
    trend = []; season = [];
    y_deseason = y;
    fitresult_linear = []; fitresult_season = [];
end

% Grid search
% correl = 0;
% init = grid_search(lb, ub, n_lag, ncontracts, LT, correl);
% init_idx = find(init(:,1) < init(:,4));
% init(init_idx,:) = [];
% log_L = kf_v1(init(1,:), par_names, y_deseason, deltat, ttm, LT, correl, serial);
% % parpool(8);
% parfor i = 1:size(init,1)
%     log_L = kf_v1(init(i,:), par_names, y_deseason, deltat, ttm, LT, correl, serial);
%     logL(i) = log_L;
% end
% par_init = mean(init(find(logL == min(logL)),:),1);
% 
% for i = 1:length(par_init)
%     if lb(i) >= par_init(i)
%         par_init(i) = par_init(i) + 1e-4;
%     elseif par_init(i) >= ub(i)
%         par_init(i) = par_init(i) - 1e-4;
%     end
% end
% 
% % if correlation == 1
% rho_range = linspace(-0.9999, 0.9999, 5);
% [grids{1:ncontracts}] = ndgrid(rho_range);
% parameter_grid = cell2mat(cellfun(@(x) x(:), grids, 'UniformOutput', false));
% init = [repmat(par_init,length(parameter_grid), 1) parameter_grid];
% parfor i = 1:length(init)
%     i;
%     log_L = kf_v1(init(i,:), par_names, y_deseason, deltat, ttm, LT, correlation, serial);
%     logL(i) = log_L;
% end
% % par_init = mean(init(find(logL == min(logL)),:),1);
% par_init = init(find(logL == min(logL)),:);
% par_init = par_init(1,:);

global save_ett save_ytt save_att save_vtt save_att_1 save_Ptt_1 save_Ptt save_vtt_star
% Setting options for optimisation
options = optimset('Display', 'iter', 'TolFun', 1e-5, 'TolX', 1e-5, 'MaxIter', 500, 'MaxFunEvals', 50000, 'OutputFcn', @check_complexity);
% options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'EnableFeasibilityMode', true,  'SubproblemAlgorithm', 'cg', 'Display', 'iter', 'TolFun', 1e-5, 'TolX', 1e-5, 'MaxIterations', 500, 'MaxFunctionEvaluations', 50000, 'OutputFcn', @check_complexity);
    
non_ar_param = model_par+ncontracts+ncontracts*(correlation == 1);
try
    phi = par_init(non_ar_param+1:end);
catch
    phi = [];
end

for j = 1:ncontracts
    phi_temp = phi(j:ncontracts:end);
    non_zero_indices = find(phi_temp ~= 0);
    if isempty(non_zero_indices) 
        orders(j) = 0; % No non-zero coefficients (edge case) 
    else 
        orders(j) = max(non_zero_indices);
    end
end
nonlcon = @(phi) stationarity_constraint_VAR_v3(phi, ncontracts, orders);

% Adjust bounds based on AR orders
for i = 1:ncontracts 
    for j = 1:n_lag
        idx = non_ar_param + (j - 1) * ncontracts + i; % Index in params vector 
        if orders(i) == 1 
            if j == 1 
                % Set bounds for AR(1) coefficients 
                lb(idx) = -0.9999; ub(idx) = 0.9999; 
            else 
                % Set bounds for higher-order coefficients to 0 (not used) 
                lb(idx) = 0; ub(idx) = 0; 
            end 
        end 
    end 
end

for i = 1:length(par_init)
    if lb(i) >= par_init(i)
        if par_init(i) == 0
            par_init(i) = 0;
        else
            par_init(i) = par_init(i) + 1e-4;
        end
    elseif par_init(i) >= ub(i)
        par_init(i) = par_init(i) - 1e-4;
    end
end

if LT == "GBM"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v3_arp, par_init, [], [], [], [], lb, ub, [], options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
elseif LT == "OU"
    % [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v3_arp, par_init, A, b, [], [], lb, ub, nonlcon, options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@(par) kf_v3_arp(par, par_names, y_deseason, deltat, ttm, LT, correlation, serial), par_init, A, b, [], [], lb, ub, nonlcon, options);

end

par_optim = par_optimised;
log_L_optim = log_L;
att = save_att;
ytt = save_ytt;
ett = save_ett;
vtt = save_vtt;
att_1 = save_att_1;
Ptt = save_Ptt;
Ptt_1 = save_Ptt_1;

% ----------------------------------------------------------
function stop = check_complexity(x, optimValues, state) 
% Check if the objective function is complex 
if isnan(optimValues.fval)
    stop = true; % Stop the optimization
else
    stop = false; % Continue the optimization
end
