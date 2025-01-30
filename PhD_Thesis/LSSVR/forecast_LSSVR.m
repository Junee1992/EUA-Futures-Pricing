function output = forecast_LSSVR(y, horiz);

data = y;
[nobsn, ncontracts] = size(data);
n_temp = nobsn - horiz;
y_data_temp = data(1:n_temp,:);
for r = 1:20
    r
    clear resid y_pred y_data
    for j = 1:ncontracts
        y_data = y_data_temp(1:n_temp+r-1,j);%-season(1:end-horiz)-trend(1:end-horiz);

        for i = 1:3
            for k = 1:3
                AR_model{i,k} = arima('Constant',NaN,'ARLags',1:i,'D',0,'MALags',[], 'Distribution','Gaussian');
                MA_model{i,k} = arima('Constant',NaN,'ARLags',[],'D',0,'MALags',1:k, 'Distribution','Gaussian');
                AR_y{i,k} = estimate(AR_model{i,k}, y_data,'Display','off');
                MA_y{i,k} = estimate(MA_model{i,k}, y_data,'Display','off');
                AIC_AR(i,k) = summarize(AR_y{i,k}).AIC;
                AIC_MA(i,k) = summarize(MA_y{i,k}).AIC;
                ARMA_model{i,k} = arima('Constant',NaN,'ARLags',1:i,'D',0,'MALags',1:k, 'Distribution','Gaussian');
                try
                    ARMA_y{i, k} = estimate(ARMA_model{i, k}, y_data, 'Display', 'off');
                    AIC_ARMA(i,k) = summarize(ARMA_y{i,k}).AIC;
                catch
                    warning('Estimation failed for ARMA_model{%d, %d}. Skipping...', i, k);
                    ARMA_y{i, k} = NaN; % or set to an alternative default value
                    AIC_ARMA(i,k) = NaN;
                end
                % ARMA_y{i,k} = estimate(ARMA_model{i,k}, y_data,'Display','off');
                ARIMA_model{i,k} = arima('Constant',NaN,'ARLags',1:i,'D', 1, 'MALags',1:k,'Distribution','Gaussian');
                try
                    ARIMA_y{i,k} = estimate(ARIMA_model{i,k},y_data,'Display','off');
                    AIC_ARIMA(i,k) = summarize(ARIMA_y{i,k}).AIC;
                catch
                    warning('Estimation failed for ARIMA_model{%d,1,%d}. Skipping...', i, k);
                    ARMA_y{i, k} = NaN; % or set to an alternative default value
                    AIC_ARIMA(i,k) = NaN;
                end
                % AIC_ARIMA(i,k) = summarize(ARIMA_y{i,k}).AIC;
            end
        end
        % info_crit = [AIC_AR; AIC_MA; AIC_ARMA];
        info_crit = [AIC_AR; AIC_MA; AIC_ARMA; AIC_ARIMA];
        info_crit(1:9,:) = NaN;
        [~,ind] = min(info_crit(:));
        [row,col] = ind2sub(size(info_crit), ind);
        if row < 4
            k = row; i = row;
            mdl{r}{j} = AR_y{i,k};
        elseif row < 7 & row > 3
            k = col; i = row-3;
            mdl{r}{j} = MA_y{i,k};
        elseif row < 10 & row > 6
            i = row-6; k = col;
            mdl{r}{j} = ARMA_y{i,k};
        elseif row > 9
            i = row-9; k = col;
            mdl{r}{j} = ARIMA_y{i,k};
        end

        resid(:,j) = infer(mdl{r}{j}, y_data);
    end

    % ncontracts = size(resid,2);
    n_lag = 3; y = resid;
    x_train = windowize(resid(1:end-1,:), 1:n_lag);
    y_train = resid(1+n_lag:end,:);
    % res = resid(4:end,:);

    % trainingX = [x_train res];
    trainingX = x_train;
    [gamma, p, MSE] = GridMLSSVR(trainingX, y_train, 5);
    % gamma = 10000; p = 10000;
    [alpha, b] = MLSSVRTrain(trainingX, y_train, gamma, p);
    % predictY = MLSSVRPredict(trainingX, y_train, trainingX, alpha, b, p);
    % plot(predictY(:,1))
    testY = zeros(horiz,ncontracts);
    fit_resid = MLSSVRPredict(trainingX, y_train, trainingX, alpha, b, p);
    res_lssvr_var(r,:) = var(y_train - fit_resid);
    testingX = [trainingX(end, 1+ncontracts:n_lag*ncontracts), y_train(end,:)];
    [e_predict] = MLSSVRPredict(testingX, testY(i,:), trainingX, alpha, b, p);
    for j = 1:ncontracts
        y_fore(:,j) = forecast(mdl{r}{j}, 1, 'Y0', y_data_temp(:,j)) + e_predict(j);
    end

    y_data_temp = [y_data_temp; y_fore];
    plot(data(:,1), 'k');
    hold on
    plot(y_data_temp(:,1), 'r');
    pause(0.2);
end
% name = sprintf('lssvm_forecast_period_%d_newnew', q);
% save(name, 'mdl', 'd ata', 'y_data_temp', 'y_fore', 'res_lssvr_var');
output = struct('y_data_temp', y_data_temp, 'res_lssvr_var', res_lssvr_var);
end
% name = sprintf('lssvm_forecast_update_%d', p);
% save('lssvm_forecast_update_v6.mat', 'my3', 'lb3', 'ub3', 'y_forecasting', 'model_fit')