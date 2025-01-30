function [output] = forecast_NPF(y, ttm, model_options)

fields = fieldnames(model_options);
for i = 1:numel(fields)
    fieldName = fields{i};
    eval([fieldName ' = model_options.(fieldName);']);
end

if LT == "OU"
    model_par = 8;
elseif LT == "GBM"
    model_par = 7;
end

n_par = length(par_names);
s2M = [repelem(2, model_par), repelem(1, n_par-model_par)]/M^(1.5); % Variances for jittering kernel
% Appendix C of Crisan and Miguez (2018)

if err == "laplace"
    n_beta = 1;
    s2M = [s2M 2/M^(1.5)];
elseif err == "hyperbolic"
    n_beta = 3+ncontracts;
    s2M = [s2M repelem(2, 3)/M^(1.5) repelem(1, ncontracts)/M^(1.5)];
end

[~, ~, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, 0);
ub(9:8+ncontracts) = 10;
lb(5) = -3; ub(5) = 10;

% Initialisation for the parameters (external particle)
n_temp = nobsn - n_forecast;
[par_optim, log_L_optim, par_init, trend, season, att, ytt, ett, ett_1, Ptt, fitresult_linear, fitresult_season] = ...
    par_estimate(y(1:n_temp,:), ttm(1:n_temp,:), model_options);
if detrend_price == "yes"
    y_deseason = y(1:nobsn-n_forecast,:) - trend - season;
else
    y_deseason = y;
end

par_optimised = par_optim;
model_options_a = struct('LT', LT, 'correlation', correlation, 'par_names', par_names,...
    'M', M, 'N', N, 'deltat', deltat, 'err', err);

par_particle_temp = [par_optim*0.7; par_optim*1.3]';
par_particle = par_particle_temp(:,1)' + (par_particle_temp(:,2)-par_particle_temp(:,1))'.*rand(M,size(par_names,2)-n_beta);
for i = 1:length(par_particle(:,1))
    while any(par_particle(i,:) > ub) || any(par_particle(i,:) < lb)
        par_particle(i,:) = par_particle_temp(:,1)' + (par_particle_temp(:,2)-par_particle_temp(:,1))'.*rand(1,size(par_particle_temp,1));
    end
end

if correlation == 1
    for i = 1:M
        % for j = 9+ncontracts:8+2*ncontracts
        for j = 9+ncontracts:length(par_names)-n_beta
            if par_particle(i,j) > 1
                par_particle(i,j) = par_particle_temp(j,1)' + (1-par_particle_temp(j,1))'.*rand(1,1);
            end
        end
    end
end

err = model_options_a.err;
if err == "laplace"
    par_particle(:,8+(1+(correlation == 1)) * ncontracts+1) = repelem(0, M);
    theta_initial = 3;

    % Lower and upper bounds for theta
    lb_s = 0.001; ub_s = 10;
    if correlation == 0
        cov_matrix = diag(par_optim(9:8+ncontracts).^2);
    else
        correl = par_optim(9+ncontracts:8+2*ncontracts);

        % Manually creating the correlation matrix of measurement errors
        CorMat = diag(repelem(1, ncontracts));
        for i = 1:ncontracts
            for j = 1:ncontracts
                if i == j
                    CorMat(i,j) = 1;
                else
                    CorMat(i,j) = correl(i) * correl(j);
                end
            end
        end
        D = diag(par_optim(9:8+ncontracts).^2);
        cov_matrix = D^(1/2) * CorMat * D^(1/2);
    end

    % Optimize both alpha and diagonal elements
    options = optimoptions('fmincon', 'Display', 'iter');
    theta_opt = fmincon(@(theta) objectiveFunction(theta, ett, cov_matrix), theta_initial, [], [], [], [], lb_s, ub_s, [], options);
    theta_opt = 2;
    theta_temp = [theta_opt * 0.5; theta_opt * 1.5]';
    par_particle(:,end) = theta_temp(:,1)' + (theta_temp(:,2)-theta_temp(:,1))'.*rand(M,size(theta_opt,2));
    if correlation == 1
        for i = 1:M
            for j = 9+ncontracts:8+2*ncontracts
                if par_particle(i,j) > 1
                    par_particle(i,j) = theta_temp(j-8,1)' + (1-theta_temp(j-8,1))'.*rand(1,1);
                end
            end
        end
    end
    lb = [lb 1e-5]; ub = [ub 3];
elseif err == "hyperbolic"
    for i = 1:ncontracts
        gmskew(i) = sprintf("nu_%d", i);
    end
    par_names = [par_names "lambda", "psi", "chi", gmskew];
    par_particle(:,8+(1+(correlation == 1)) * ncontracts+1:8+(1+(correlation == 1)) * ncontracts+ncontracts+3) = zeros(M, 3+ncontracts);
    % theta_initial = [(par_optim(9:8+(1+(correlation == 1)) * ncontracts)) 3];  % Initial guesses for alpha and the standard deviations
    % Lower and upper bounds for theta
    % lb_s = [0.00001*ones((1+(correlation == 1)) * ncontracts,1); 0.001];  % Lower bound for alpha and the standard deviations (assumed bounds)
    % ub_s = [1*ones((1+(correlation == 1)) * ncontracts,1); 10];  % Upper bound for alpha and the standard deviations
    theta_initial = [1.5, 1.5, 1.5, repelem(0.005, ncontracts)];
    lb_s = [repelem(0.001, 3), repelem(-1e-2, ncontracts)]; ub_s = [repelem(10, 3), repelem(1e-2, ncontracts)];
    if correlation == 0
        cov_matrix = diag(par_optim(9:8+ncontracts).^2);
    else
        correl = par_optim(9+ncontracts:8+2*ncontracts);

        % Manually creating the correlation matrix of measurement errors
        CorMat = diag(repelem(1, ncontracts));
        for i = 1:ncontracts
            for j = 1:ncontracts
                if i == j
                    CorMat(i,j) = 1;
                else
                    CorMat(i,j) = correl(i) * correl(j);
                end
            end
        end
        D = diag(par_optim(9:8+ncontracts).^2);
        cov_matrix = D^(1/2) * CorMat * D^(1/2);
    end
    % Optimize both alpha and diagonal elements
    options = optimoptions('fmincon', 'Display', 'iter');
    theta_opt = fmincon(@(theta) GHobjectiveFunction(theta, ett, cov_matrix), theta_initial, [], [], [], [], lb_s, ub_s, [], options);
    theta_temp = [theta_opt' * 0.7, theta_opt' * 1.3];
    % par_particle(:,8+2*ncontracts+1) = gamrnd(1, 100, M, 1);
    par_particle(:,8+(1+(correlation == 1)) * ncontracts+1:end) = theta_temp(:,1)' + (theta_temp(:,2)-theta_temp(:,1))'.*rand(M,size(theta_opt,2));
    if correlation == 1
        for i = 1:M
            for j = 9+ncontracts:8+2*ncontracts
                if par_particle(i,j) > 1
                    par_particle(i,j) = theta_temp(j-8,1)' + (1-theta_temp(j-8,1))'.*rand(1,1);
                end
            end
        end
    end
    % if err == "laplace"
    % lb = [lb 1e-5]; ub = [ub 3];
    % elseif err == "hyperbolic"
    % lb = [lb -1 1e-5 1e-5 repelem(-1e-2, ncontracts)];
    lb = [lb 1e-3 1e-3 1e-3 repelem(-1e-2, ncontracts)];
    ub = [ub 10 10 10 repelem(1e-2, ncontracts)];
    % end
end
% y = ett;
%
for i = 1:length(par_particle(:,1))
    while par_particle(i,1) < par_particle(i,4)
        % par_particle(i, 4) = betarnd(1.5, 5); % kappa >= gamma, Binkowski et al (2019)
        par_particle(i,4) = par_particle_temp(4,1)' + (par_particle_temp(4,2)-par_particle_temp(4,1))'.*rand(1,1);
    end

end
part_est = zeros(nobsn-n_forecast, 2); % estimates of states
par_est = zeros(nobsn-n_forecast, size(par_names,2)); % estimates of parameters
% y = y(:,[1 3 4]);

%
% [~, att_up, ~, ytt_up, ~] = kf_v1_laplace([par_optim(1:8), 100*par_optim(9:15), theta_opt], ett, ttm, model_options);
% part_est(1,:) = [0.0001, x(1,2)];
particle = cell(N, 1);
for i = 1:M
    particle{i} = [att(1,1)*ones(N,1), att(1,2) * ones(N, 1)] + 0.01 * randn(N, 2);
    % particle{i} = [att_up(1,1)*ones(N,1), att_up(1,2) * ones(N, 1)] + 0.001 * randn(N, 2);
    % particle{i} = [x(1,1)*ones(N,1), x(1,2)*ones(N,1)] + 0.1 * randn(N,2);
end

lq_chi = zeros(nobsn-n_forecast, 1);
uq_chi = lq_chi;
uq_xi = lq_chi;
lq_xi = lq_chi;
part_post = zeros(M, size(particle{i}, 2));
% rng(10000);
% y_temp = y(1:100,:) - season-trend;
y_temp = y_deseason;
part_est(1,:) = [att(1,1), att(1,2)];
% part_est(1,:) = [att_up(1,1), att_up(1,2)];

%
for t = 1:nobsn-n_forecast
    %for t= 1:1
    y_t = y_temp(t,:); ttm_t = ttm(t,:);
    lW = zeros(1, M);
    particle_post_cell = cell(100,1);
    % figure;
    % yline(y(1,1));
    % Parallelise
    parfor i = 1:M
        [t, i]
        if t > 1
            local_particle = par_particle(i, :); % Create a local copy to avoid conflicts
            for q = 1:length(local_particle)
                if q == 4 % setting kappa >= gamma.
                    if local_particle(1) < local_particle(q)
                        local_particle(q) = rnd_tgaussian(local_particle(q), s2M(q)*ones(1), lb(q), local_particle(1));
                    else
                        local_particle(q) = rnd_tgaussian(local_particle(q), s2M(q)*ones(1), lb(q), ub(q));
                    end
                else
                    local_particle(q) = rnd_tgaussian(local_particle(q), s2M(q)*ones(1), lb(q), ub(q));
                end
            end
            par_particle(i, :) = local_particle;
        end

        % Inner filter
        [lW(i), particle_post_cell{i}, particle{i}, wu{i}, w{i}, n_eff]  = kf_onestep_neff_v2(y_t, ttm_t, par_particle(i,:), particle{i}, model_options_a, part_est, i);
        part_post(i,:) = particle_post_cell{i};
        wu_save(i,:) = wu{i};
        w_save(i,:) = w{i};
        n_eff_int(i) = n_eff;
    end
    %
    % subplot(2,2,4)
    % plot(n_eff_int, '*');
    % hold on; pause(0.01)
    n_eff_int_total{t} = n_eff_int;
    % n_eff_int_total(t,:) = n_eff_int;
    % weight_unnormalised{t} = wu_save;
    % weight_normalised{t} = w_save;

    % Calculate weights for external particles
    w_o = exp(lW-max(lW));
    w_norm = w_o/sum(w_o);
    part_est(t,:) = (part_post'* w_norm')';
    par_est(t,:) = (par_particle' * w_norm')';
    ext_weight_unnormalised{t} = w_o;
    ext_weight_normalised{t} = w_norm;
    % Resampling
    Neff(t) = 1/sum(w_norm.^2);
    if Neff(t)<0.3*M
        [mu, si, p] = EMGM(par_particle',w_norm',5);
        nw=length(p);
        % nw_i=[nw_i,nw];
        for j=1:M
            indw  = randsample(nw, 1, true,p);
            par_particle(j,:) = mvnrnd(mu(:,indw), si(:,:,indw));
            while any(par_particle(j,:) > ub) || any(par_particle(j,:) < lb)
                par_particle(j,:) = mvnrnd(mu(:,indw), si(:,:,indw));
            end
        end
        % w=1/M*ones(M,1); % reweighting
        % H2=-w'*log(mvnpdf(part_u,mu_u,cov_u)+realmin);
        % H=[H;H1,H2];

    end

    % if Neff(t)<0.5*M
    % [par_particle_resamp, ~] = Resample(par_particle', w_norm);
    % par_particle = par_particle_resamp';
    % end

    % Intervals
    % Short-term factor
    data = [part_post(:,1), w_norm'];
    sorted_data = sortrows(data, 1);
    cumprob = cumsum(sorted_data(:,2));

    [~,ind] = min(abs(cumprob-0.025));
    lq_chi(t,:) = sorted_data(ind,1);

    [~,ind] = min(abs(cumprob-0.975));
    uq_chi(t,:) = sorted_data(ind,1);

    % Long-term factor
    data = [part_post(:,2), w_norm'];
    sorted_data = sortrows(data, 1);
    cumprob = cumsum(sorted_data(:,2));

    [~,ind] = min(abs(cumprob-0.025));
    lq_xi(t,:) = sorted_data(ind,1);

    [~,ind] = min(abs(cumprob-0.975));
    uq_xi(t,:) = sorted_data(ind,1);

    subplot(1,2,1)
    plot(1:t, Neff(1:t)/M);
    subplot(1,2,2)
    plot(1:t, par_est(1:t,:));
    pause(0.01);
end

bound_chi = [lq_chi uq_chi];
bound_xi = [lq_xi uq_xi];
param_est = par_est;
state_est = part_est;

% Forecast
serial = "no"; n_lag = 0;
model_options_b = struct('LT', LT, 'correlation', correlation, 'par_names', par_names, 'M', M, 'N', N, 'deltat', deltat, 'err', err, 'Serialcor', serial);
par = set_parameters(LT, ncontracts, par_names, param_est(end,:), correlation, "no", err);

C_test = [0; par.mu * (1 - exp(-par.gamma * deltat)) / par.gamma];
G_test = [exp(-par.kappa * deltat) , 0; 0, exp(-par.gamma * deltat)];
% x_test = zeros(size(C_test,1), round(n_forecast)+1);
% if p == 1
%     a0 = att(end-20,:)';
% elseif p == 2
% att = att_set{p};
% if serial == "no"
%     att = att;
% end
a0 = state_est(end,:)';
a0_save(:,1) = a0;
x_test(1,:) = C_test + G_test * a0;
[~, ~, Ptt] = kf_v5_nonG(param_est(end,:), y_deseason(1:n_temp,:), ttm(1:n_temp,:), model_options_b);

Ptt_save(:,:,1) = Ptt(:,:,end);
% end
% x_test(:,1) = att(end,:)';

W = [(1 - exp(-2 * par.kappa * deltat)) / (2 * par.kappa) * par.sigmachi^2, (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi);
    (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi), (1 - exp(-2 * par.gamma * deltat)) / (2 * par.gamma) * par.sigmaxi^2];
ttm_test = ttm(n_temp-n_lag+1:n_temp+1, :);

% Measurement Equation
d1_test = (1 - exp(-2 * par.kappa * ttm_test)) * par.sigmachi^2 / (2 * par.kappa);
d2_test = (1 - exp(-2 * par.gamma * ttm_test)) * par.sigmaxi^2 / (2 * par.gamma);
d3_test = (1 - exp(-(par.kappa + par.gamma) * ttm_test)) * 2 * par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa + par.gamma);
d_temp = (par.mu - par.lambdaxi) / par.gamma * (1 - exp(-par.gamma * ttm_test)) -...
    (par.lambdachi / par.kappa) * (1 - exp(-par.kappa * ttm_test)) + ...
    (1/2) * (d1_test + d2_test + d3_test);
d_temp = d_temp';

for i = 1:size(ttm_test,1)
    B1(i,:) = exp(-par.kappa * ttm_test(i,:));
    B2(i,:) = exp(-par.gamma * ttm_test(i,:));
    B_temp(:,:,i) = [B1(i,:); B2(i,:)]';
end

if correlation == 0
    V_temp = diag(par.s.^2);

elseif correlation == 1
    correl = par.rho;

    % Manually creating the correlation matrix of measurement errors4
    CorMat = diag(repelem(1, ncontracts));
    for i = 1:ncontracts
        for j = 1:ncontracts
            if i == j
                CorMat(i,j) = 1;
            else
                CorMat(i,j) = correl(i) * correl(j);
            end
        end
    end
    D = diag(par.s.^2);
    V_temp = D^(1/2) * CorMat * D^(1/2);

else
    error('correlation must be 0 or 1.')
end
V_temp = chol(V_temp)'*chol(V_temp);
% fitresult_season = season_set{p}; fitresult_linear = trend_set{p};

y_temp = y(1:n_temp,:);

% if serial == "no"
    d_test = d_temp; B_test = B_temp;
    if err == "laplace"
        % for j = 1:ncontracts
        %     v_test(j) = gal_inv_cdf(0.975, V_temp(j,j), par.sG);
        % end
        y_test(1,:) = d_test + B_test * x_test(1,:)';
        varn(:,:,1) = B_test * (G_test * Ptt(:,:,end) * G_test' + W) * B_test' + par.sG.*V_temp;
    elseif err == "hyperbolic"
        for i = 1:10000
            uv(i,:) = gigrnd(par.lambda, par.psi, par.chi).*par.nu;
        end
        y_test(1,:) = d_test + B_test * x_test(1,:)' + mean(uv)';
        var_gig = (par.chi/par.psi)*besselk(par.lambda+2, sqrt(par.chi*par.psi))/besselk(par.lambda, sqrt(par.chi*par.psi)) - ...
            ((par.chi/par.psi).^0.5 * besselk(par.lambda+1, sqrt(par.chi*par.psi))/besselk(par.lambda, sqrt(par.chi*par.psi))).^2;
        e_gig = (par.chi/par.psi).^0.5 * besselk(par.lambda+1, sqrt(par.chi*par.psi))/besselk(par.lambda, sqrt(par.chi*par.psi));
        var_gh(:,:,1) = par.nu'*par.nu*var_gig + e_gig*V_temp;
        varn(:,:,1) = B_test * (G_test * Ptt(:,:,end) * G_test' + W) * B_test' + var_gh(:,:,1);
    end
    if detrend_price == "yes"
        y_fore(1,:) = y_test(1,:) + fitresult_season(n_temp+1) + fitresult_linear(n_temp+1);
    else
        y_fore(1,:) = y_test(1,:);
    end
    y_temp = [y_temp; y_fore(1,:)];
   
% end

for t = n_temp+1:n_temp+n_forecast-1
    %for t= 1:1
    if detrend_price == "yes"
        % Detrend the price
        [y_detrend, trend, fitresult_linear] = linear_detrend(y_temp);
        % [season, y_deseason, fitresult_season, gof] = deseasonalise_v2(y_detrend);
        seasonality = "sinusoid_b";
        for i = 1:10
            [season, p, fitresult_season, gof] = deseasonalise(y_detrend, i, seasonality);
            rmse_seasonality(i) = gof.rmse;
        end
        no_season = min(find(rmse_seasonality == min(rmse_seasonality)));
        [season, p, fitresult_season, gof] = deseasonalise(y_detrend, no_season, seasonality);
        y_d= y_detrend - season;
    else
        trend = []; season = []; fitresult_linear = []; fitresult_season = [];
        y_d= y;
    end

    y_t = y_d(t,:); ttm_t = ttm(t,:);
    lW = zeros(1, M);
    particle_post_cell = cell(100,1);
    % figure;
    % yline(y(1,1));
    % Parallelise
    parfor i = 1:M
        [t, i]
        % if t > 1
        %     local_particle = par_particle(i, :); % Create a local copy to avoid conflicts
        %     for q = 1:length(local_particle)
        %         if q == 4 % setting kappa >= gamma.
        %             if local_particle(1) < local_particle(q)
        %                 local_particle(q) = rnd_tgaussian(local_particle(q), s2M(q)*ones(1), lb(q), local_particle(1));
        %             else
        %                 local_particle(q) = rnd_tgaussian(local_particle(q), s2M(q)*ones(1), lb(q), ub(q));
        %             end
        %         else
        %             local_particle(q) = rnd_tgaussian(local_particle(q), s2M(q)*ones(1), lb(q), ub(q));
        %         end
        %     end
        %     par_particle(i, :) = local_particle;
        % end

        % Inner filter
        [lW(i), particle_post_cell{i}, particle{i}, wu{i}, w{i}, n_eff]  = kf_onestep_neff_v2(y_t, ttm_t, par_particle(i,:), particle{i}, model_options_b, part_est, i);
        part_post(i,:) = particle_post_cell{i};
        wu_save(i,:) = wu{i};
        w_save(i,:) = w{i};
        n_eff_int(i) = n_eff;
    end
    %
    % subplot(2,2,4)
    % plot(n_eff_int, '*');
    % hold on; pause(0.01)
    n_eff_int_total{t} = n_eff_int;
    % n_eff_int_total(t,:) = n_eff_int;
    % weight_unnormalised{t} = wu_save;
    % weight_normalised{t} = w_save;

    % Calculate weights for external particles
    w_o = exp(lW-max(lW));
    w_norm = w_o/sum(w_o);
    part_est(t,:) = (part_post'* w_norm')';
    par_est(t,:) = (par_particle' * w_norm')';
    ext_weight_unnormalised{t} = w_o;
    ext_weight_normalised{t} = w_norm;
    % Resampling
    Neff(t) = 1/sum(w_norm.^2);
    if Neff(t)<0.3*M
        [mu, si, p] = EMGM(par_particle',w_norm',5);
        nw=length(p);
        % nw_i=[nw_i,nw];
        for j=1:M
            indw  = randsample(nw, 1, true,p);
            par_particle(j,:) = mvnrnd(mu(:,indw), si(:,:,indw));
            while any(par_particle(j,:) > ub) || any(par_particle(j,:) < lb)
                par_particle(j,:) = mvnrnd(mu(:,indw), si(:,:,indw));
            end
        end
        % w=1/M*ones(M,1); % reweighting
        % H2=-w'*log(mvnpdf(part_u,mu_u,cov_u)+realmin);
        % H=[H;H1,H2];

    end

    % if Neff(t)<0.5*M
    % [par_particle_resamp, ~] = Resample(par_particle', w_norm);
    % par_particle = par_particle_resamp';
    % end

    % Intervals
    % Short-term factor
    data = [part_post(:,1), w_norm'];
    sorted_data = sortrows(data, 1);
    cumprob = cumsum(sorted_data(:,2));

    [~,ind] = min(abs(cumprob-0.025));
    lq_chi(t,:) = sorted_data(ind,1);

    [~,ind] = min(abs(cumprob-0.975));
    uq_chi(t,:) = sorted_data(ind,1);

    % Long-term factor
    data = [part_post(:,2), w_norm'];
    sorted_data = sortrows(data, 1);
    cumprob = cumsum(sorted_data(:,2));

    [~,ind] = min(abs(cumprob-0.025));
    lq_xi(t,:) = sorted_data(ind,1);

    [~,ind] = min(abs(cumprob-0.975));
    uq_xi(t,:) = sorted_data(ind,1);
    param_est = par_est;
    state_est = part_est;

    par = set_parameters(LT, ncontracts, par_names, param_est(end,:), correlation, "no", err);

    C_test = [0; par.mu * (1 - exp(-par.gamma * deltat)) / par.gamma];
    G_test = [exp(-par.kappa * deltat) , 0; 0, exp(-par.gamma * deltat)];
    % x_test = zeros(size(C_test,1), round(n_forecast)+1);
    % if p == 1
    %     a0 = att(end-20,:)';
    % elseif p == 2
    % att = att_set{p};
    % if serial == "no"
    %     att = att;
    % end
    a0 = state_est(end,:)';
    a0_save(:,t-n_temp+1) = a0;
    x_test(t-n_temp+1,:) = C_test + G_test * a0;
    [~, ~, Ptt] = kf_v5_nonG(param_est(end,:), y_d(1:t,:), ttm(1:t,:), model_options_b);
    Ptt_save(:,:,t-n_temp+1) = Ptt(:,:,end);
    % end
    % x_test(:,1) = att(end,:)';

    W = [(1 - exp(-2 * par.kappa * deltat)) / (2 * par.kappa) * par.sigmachi^2, (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi);
        (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi), (1 - exp(-2 * par.gamma * deltat)) / (2 * par.gamma) * par.sigmaxi^2];
    ttm_test = ttm(t+1, :);

    % Measurement Equation
    d1_test = (1 - exp(-2 * par.kappa * ttm_test)) * par.sigmachi^2 / (2 * par.kappa);
    d2_test = (1 - exp(-2 * par.gamma * ttm_test)) * par.sigmaxi^2 / (2 * par.gamma);
    d3_test = (1 - exp(-(par.kappa + par.gamma) * ttm_test)) * 2 * par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa + par.gamma);
    d_temp = (par.mu - par.lambdaxi) / par.gamma * (1 - exp(-par.gamma * ttm_test)) -...
        (par.lambdachi / par.kappa) * (1 - exp(-par.kappa * ttm_test)) + ...
        (1/2) * (d1_test + d2_test + d3_test);
    d_temp = d_temp';

    for i = 1:size(ttm_test,1)
        B1(i,:) = exp(-par.kappa * ttm_test(i,:));
        B2(i,:) = exp(-par.gamma * ttm_test(i,:));
        B_temp(:,:,i) = [B1(i,:); B2(i,:)]';
    end

    if correlation == 0
        V_temp = diag(par.s.^2);

    elseif correlation == 1
        correl = par.rho;

        % Manually creating the correlation matrix of measurement errors4
        CorMat = diag(repelem(1, ncontracts));
        for i = 1:ncontracts
            for j = 1:ncontracts
                if i == j
                    CorMat(i,j) = 1;
                else
                    CorMat(i,j) = correl(i) * correl(j);
                end
            end
        end
        D = diag(par.s.^2);
        V_temp = D^(1/2) * CorMat * D^(1/2);

    else
        error('correlation must be 0 or 1.')
    end
    V_temp = chol(V_temp)'*chol(V_temp);
    % fitresult_season = season_set{p}; fitresult_linear = trend_set{p};

    if serial == "no"
        d_test = d_temp; B_test = B_temp;
        
        if err == "laplace"
            y_test(t-n_temp+1,:) = d_test + B_test * x_test(t-n_temp+1,:)';
            varn(:,:,t-n_temp+1) = B_test * (G_test * Ptt(:,:,end) * G_test' + W) * B_test' + par.sG.*V_temp;
        elseif err == "hyperbolic"
            for i = 1:10000
                uv(i,:) = gigrnd(par.lambda, par.psi, par.chi).*par.nu;
            end
            y_test(t-n_temp+1,:) = d_test + B_test * x_test(t-n_temp+1,:)' + mean(uv)';
            var_gig = (par.chi/par.psi)*besselk(par.lambda+2, sqrt(par.chi*par.psi))/besselk(par.lambda, sqrt(par.chi*par.psi)) - ...
                ((par.chi/par.psi).^0.5 * besselk(par.lambda+1, sqrt(par.chi*par.psi))/besselk(par.lambda, sqrt(par.chi*par.psi))).^2;
            e_gig = (par.chi/par.psi).^0.5 * besselk(par.lambda+1, sqrt(par.chi*par.psi))/besselk(par.lambda, sqrt(par.chi*par.psi));
            var_gh(:,:,t-n_temp+1) = par.nu'*par.nu*var_gig + e_gig*V_temp;
            varn(:,:,t-n_temp+1) = B_test * (G_test * Ptt(:,:,end) * G_test' + W) * B_test' + var_gh(:,:,r);
        end
        
        if detrend_price == "yes"
            y_fore(t-n_temp+1,:) = y_test(t-n_temp+1,:) + fitresult_season(t+1) + fitresult_linear(t+1);
        elseif detrend_price == "no"
            y_fore(t-n_temp+1,:) = y_test(t-n_temp+1,:);
        end

        y_temp = [y_temp; y_fore(t-n_temp+1,:)];
       

    end
    % plot(400:t+1, y(400:t+1,1), 'k')
    % hold on
    % plot(460:t+1, y_temp(460:t+1,1), 'g')
    % pause(0.2)

end
output = struct('y_temp', y_temp, 'varn', varn, 'param_est', param_est, 'state_est', state_est);
