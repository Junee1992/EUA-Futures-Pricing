function [log_L, att, Ptt, ytt, a_s] = kf_v5_nonG(par_init, y, ttm, model_options)

fields = fieldnames(model_options);
for i = 1:numel(fields)
    fieldName = fields{i};
    eval([fieldName ' = model_options.(fieldName);']);
end

serial = Serialcor;

% Defining the dimensions of the dataset.
[nobsn, ncontracts] = size(y);
n_par = size(par_init, 2);
par = set_parameters(LT, ncontracts, par_names, par_init, correlation, serial, err);
if serial == "yes"
    n_lags = length(par.phi) / ncontracts;
end

% Transition matrices in state and measurement equation
% x_t = C + G x_{t-1} + w_t,  w_t ~ N(0, W)
% y_t = d_t + B_t * x_t + v_t, v_t ~ N(0, V)
if LT == "GBM"

    C = [0; par.mu * deltat];
    G = [exp(-par.kappa * deltat), 0; 0, 1];
    W = [(1 - exp(-2 * par.kappa * deltat)) / (2 * par.kappa) * par.sigmachi^2, (1 - exp(-par.kappa * deltat)) / par.kappa * (par.sigmachi * par.sigmaxi * par.rho_chixi);...
        (1 - exp(-par.kappa * deltat)) / par.kappa * (par.sigmachi * par.sigmaxi * par.rho_chixi), par.sigmaxi^2 * deltat]; %the covariance matrix of w

    d1 = (1 - exp(-2 * par.kappa * ttm)) * par.sigmachi^2 / (2 * par.kappa);
    d2 = par.sigmaxi^2 * ttm;
    d3 = (1 - exp(-par.kappa * ttm)) * 2 * par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa);
    d = (par.mu - par.lambdaxi) * ttm - (par.lambdachi / par.kappa) * (1 - exp(-par.kappa * ttm)) + (1/2) * (d1 + d2 + d3);
    d = d';

    for i = 1:nobsn
        B1(i,:) = exp(-par.kappa * ttm(i,:));
        B2(i,:) = repelem(1, ncontracts);
        B(:,:,i) = [B1(i,:); B2(i,:)]';
    end

elseif LT == "OU"

    C = [0; (par.mu / par.gamma) * (1 - exp(-par.gamma * deltat))];
    G = [exp(-par.kappa * deltat), 0; 0, exp(-par.gamma * deltat)];
    W = [(1 - exp(-2 * par.kappa * deltat)) / (2 * par.kappa) * par.sigmachi^2, (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi);
        (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi), (1 - exp(-2 * par.gamma * deltat)) / (2 * par.gamma) * par.sigmaxi^2];

    d1 = (1 - exp(-2 * par.kappa * ttm)) * par.sigmachi^2 / (2 * par.kappa);
    d2 = (1 - exp(-2 * par.gamma * ttm)) * par.sigmaxi^2 / (2 * par.gamma);
    d3 = (1 - exp(-(par.kappa + par.gamma) * ttm)) * 2 * par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa + par.gamma);
    d = (par.mu - par.lambdaxi) / par.gamma * (1 - exp(-par.gamma * ttm)) - (par.lambdachi / par.kappa) * (1 - exp(-par.kappa * ttm)) + (1/2) * (d1 + d2 + d3);
    d = d';

    for i = 1:nobsn
        B1(i,:) = exp(-par.kappa * ttm(i,:));
        B2(i,:) = exp(-par.gamma * ttm(i,:));
        B(:,:,i) = [B1(i,:); B2(i,:)]';
    end

else
    error('Please specify the process of the long-term factor LT.')
end

vsig2 = par.s.^2;

if correlation == 0
    V = diag(vsig2);

elseif correlation == 1
    correl = par.rho;

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
    D = diag(vsig2);
    V = D^(1/2) * CorMat * D^(1/2);

else
    error('correlation must be 0 or 1.')
end

% Kalman Filter
% Initial distribution of state variables
% E(x_{1|0}) = a0
% Cov(x_{1|0}) = P0
if LT == "GBM"
    a0 = [0; mean(y(1,:), 2)];
    P0 = [0.01, 0.01; 0.01, 0.01];

elseif LT == "OU"
    a0 = [0; mean(y(1,:), 2)];
    P0 = 0.01*eye(2);
end

global save_ytt save_att save_ett save_vt save_att_1 save_Ptt_1 save_Ptt
save_ytt = zeros(nobsn, ncontracts); % Estimation of y given x_1:t, y_{t-1}
save_att = zeros(nobsn, size(a0, 1)); % Updated estimation of x at time t
save_ett = zeros(nobsn, ncontracts); % Prediction error
save_vt = zeros(nobsn, ncontracts); % Measurement errors
save_att_1 = zeros(nobsn, size(a0, 1));
save_Ptt_1 = zeros(size(a0, 1), size(a0, 1), nobsn);
save_Ptt = zeros(size(a0, 1), size(a0, 1), nobsn);

% Initial state prediction
att_1 = a0;
Ptt_1 = P0;
save_att_1(1, :) = att_1;
save_Ptt_1(:, :, 1) = Ptt_1;

negLogLikelihood = 0; % Initialize negative log-likelihood

% Importance sampling
num_samples = 100;  % Number of samples for importance sampling
% alpha_q = par.sG;  % Parameter for the proposal distribution
% beta_q = 1;     % Typically set to 1 for simplicity
if err == "laplace"
     w_samples = gamrnd(par.sG, 1, num_samples, 1);  % Draw samples from Gamma(alpha_q, beta_q)
elseif err == "hyperbolic"
    for r = 1:num_samples
       w_samples(r) = gigrnd(par.lambda, par.psi, par.chi);
    end
end

for t = 1:nobsn
    yt = y(t, :)';
   
    likelihoods = zeros(num_samples, 1);
    att_sum = zeros(size(att_1));
    Ptt_sum = zeros(size(Ptt_1));
    ytt_1 = d(:,t) + B(:,:,t) * att_1;
    
    for j = 1:num_samples
        if err == "laplace"
            et = yt - ytt_1;
        elseif err == "hyperbolic"
            et = yt - ytt_1 - w_samples(j) * par.nu';
        end
        w = w_samples(j);
        Ltt_1 = B(:, :, t) * Ptt_1 * B(:, :, t)' + w * V;

        % Regularization
        epsilon = 1e-5;
        S_j = Ltt_1 + epsilon * eye(size(Ltt_1));

        % Ensure S_j is symmetric
        S_j = (S_j + S_j') / 2;

        % Use Cholesky decomposition for numerical stability
        [L, p] = chol(S_j, 'lower');
        if p > 0
            % Further regularize if not positive definite
            S_j = S_j + 1e-5 * eye(size(S_j));
            [L, p] = chol(S_j, 'lower');
            if p > 0
                error('Covariance matrix is not positive definite even after regularization.');
            end
        end
        invS_j = L' \ (L \ eye(size(S_j)));

        % Compute innovation (prediction error) and its covariance
        % ytt_1 = d(:, t) + B(:, :, t) * att_1;
        % nu_j = yt - ytt_1;
        
        % Update step with sample w
        Kt = Ptt_1 * B(:, :, t)' * invS_j;
        att_temp(:,j) = att_1 + Kt * et;
        Ptt_temp(:,:,j) = (eye(size(Kt, 1)) - Kt * B(:, :, t)) * Ptt_1;

        % Compute likelihood
        likelihoods(j) = -0.5*(log(det(S_j)) + et' * invS_j * et);
    end

    % Normalize weights
    % weights = gampdf(w_samples, par.sG, 1) ./ gampdf(w_samples, alpha_q, beta_q);  % Compute importance weights
    weights = exp(likelihoods - max(likelihoods));
    weights = weights / sum(weights);

    % Compute weighted estimates
    for j = 1:num_samples
        att_sum = att_sum + weights(j) * att_temp(:,j);
        Ptt_sum = Ptt_sum + weights(j) * Ptt_temp(:,:,j);
    end

    % Update state estimates and covariance with weighted sums
    att = att_sum;
    Ptt = Ptt_sum;

    % Predict next state
    att_1 = C + G * att;
    Ptt_1 = G * Ptt * G' + W;

    % Compute negative log-likelihood contribution
    weighted_likelihood = sum(weights .* exp(likelihoods));
    negLogLikelihood = negLogLikelihood - log(weighted_likelihood);
    save_att(t,:) = att';
    save_ytt(t,:) = d(:,t) + B(:,:,t)*att;
    save_ett(t,:) = y(t,:) - save_ytt(t,:);

end
log_L = negLogLikelihood;

end