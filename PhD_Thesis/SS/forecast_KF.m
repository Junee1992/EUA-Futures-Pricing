function [output] = forecast_KF(y_temp, ttm, model_options)

fields = fieldnames(model_options);
for i = 1:numel(fields)
    fieldName = fields{i};
    eval([fieldName ' = model_options.(fieldName);']);
end

[~, ncontracts] = size(y_temp);

for r = 1:n_forecast
    r
    [par_optim, log_L_optim, par_init, trend, season, att, ytt, ett, ett_1, Ptt, fitresult_linear, fitresult_season] = ...
        par_estimate(y_temp(1:n_temp+r-1,:), ttm(1:n_temp+r-1,:), model_options);
    optim_param_prior{r} = par_optim;
    if LT == "GBM"
        model_par = 7;
    elseif LT == "OU"
        model_par = 8;
    end

    if length(par_optim) > model_par + 2*ncontracts
        serial = "yes";
        n_lag = (length(par_optim) - (model_par+2*ncontracts))/ncontracts;
        par_names_temp = define_parameters(LT, ncontracts, correlation, n_lag)';
        n_par_temp = length(par_names_temp);
        par_optim_temp = par_optim;
    elseif length(par_optim) == model_par + 2*ncontracts
        serial = "no";
        n_lag = 0;
        par_optim_temp = par_optim;
        par_names_temp = par_names;
        n_par_temp = length(par_names_temp);
    elseif length(par_optim) == model_par + ncontracts
        serial = "no";
        n_lag = 0;
        par_optim_temp = par_optim;
        par_names_temp = par_names;
        n_par_temp = length(par_names_temp);
    end
    if serial == "yes"
        [par_optim_temp, log_L_optim, par_init, trend, season, att, ytt, ett, vtt, fitresult_linear, fitresult_season, att_1, Ptt, Ptt_1] = ...
            param_estim_arp(y_temp(1:n_temp+r-1, :), ttm(1:n_temp+r-1, :), deltat, detrend_price, n_par_temp, par_names_temp, LT, correlation, par_optim);
    end
    optim_param{r}= par_optim_temp;

    ett_save{r} = ett; att_save{r} = att; ytt_save{r} = ytt;
    vtt_save{r} = vtt;
    lin_model{r} = fitresult_linear; seas_model{r} = fitresult_season;
    par = set_parameters(LT, ncontracts, par_names_temp, par_optim_temp, correlation, serial, "normal");

    % Parameters
    % State equation
    C_test = [0; par.mu * (1 - exp(-par.gamma * deltat)) / par.gamma];
    G_test = [exp(-par.kappa * deltat) , 0; 0, exp(-par.gamma * deltat)];
    a0 = att(end,:)';
    a0_save(:,r) = a0
    x_test(r,:) = C_test + G_test * a0;
    Ptt_save(:,:,r) = Ptt(:,:,end);

    W = [(1 - exp(-2 * par.kappa * deltat)) / (2 * par.kappa) * par.sigmachi^2, (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi);
        (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi), (1 - exp(-2 * par.gamma * deltat)) / (2 * par.gamma) * par.sigmaxi^2];
    ttm_test = ttm(n_temp-n_lag+r:n_temp+r, :);

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

    if serial == "no"
        d_test = d_temp; B_test = B_temp;
        y_test(r,:) = d_test + B_test * x_test(r,:)';
        y_fore(r,:) = y_test(r,:) + fitresult_season(n_temp+r) + fitresult_linear(n_temp+r);
        y_temp = [y_temp; y_fore(r,:)];
        varn(:,:,r) = B_test * (G_test * Ptt(:,:,end) * G_test' + W) * B_test' + V_temp;
    elseif serial == "yes"
        y_train = y_temp(1:end,:) - season - trend;
        x_temp = att(end-n_lag+1:end,:)';
        phiy = zeros(ncontracts, 1); phiB = zeros(ncontracts,1);
        phid = zeros(ncontracts, 1); phiBG = zeros(ncontracts,1);
        phiBGinv = zeros(ncontracts, size(G_test,2)); phiBGWG = zeros(ncontracts,ncontracts);
        GWG = zeros(size(G_test,1), size(G_test,2)); GBphi = zeros(size(G_test,1), ncontracts);

        for j = 1:n_lag
            phiy = phiy + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*y_train(end-j+1,:)';
            phid = phid + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*d_temp(:,n_lag+1-j);
            Ginv = 0;
            for k = 1:j
                Ginv = Ginv + inv(G_test)^k;
            end
            phiBG = phiBG + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*B_temp(:,:,end-j)*Ginv*C_test;
            phiBGinv = phiBGinv + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*B_temp(:,:,end-j)*inv(G_test)^j;
            GBphi = GBphi + (inv(G_test)')^(j-1)*B_temp(:,:,end-j)'*diag(par.phi((j-1)*ncontracts+1:j*ncontracts));

        end

        d_test = phiy + d_temp(:,end)-phid + phiBG;

        x_test_temp = C_test + G_test * a0;
        B_comb = B_temp(:,:,end) - phiBGinv;

        y_fore(r,:) = d_test(:,end) + B_comb * x_test_temp + fitresult_linear(n_temp+r) + fitresult_season(n_temp+r);
        y_temp = [y_temp; y_fore(r,:)];
        varn(:,:,r) = B_comb * Ptt_1 * B_comb'
    end
end

output = struct('y_temp', y_temp, 'y_fore', y_fore, 'optim_param_prior', optim_param_prior, 'optim_param', optim_param,...
    'ett_save', ett_save, 'ytt_save', ytt_save, 'vtt_save', vtt_save, 'att_save', att_save, 'lin_model', lin_model, ...
    'seas_model', seas_model, 'a0_save', a0_save, 'x_test', x_test, 'Ptt_save', Ptt_save, 'varn', varn);