%% Kalman filter

function [log_L] = kf_SS(par_init, y, ttm, model_options)

fields = fieldnames(model_options);
for i = 1:numel(fields)
    fieldName = fields{i};
    eval([fieldName ' = model_options.(fieldName);']);
end

% Defining the dimensions of the dataset.
[nobsn, ncontracts] = size(y);
n_par = size(par_init, 2);
par = set_parameters(LT, ncontracts, par_names, par_init, correlation, "no", "normal");

% Transition matrices in state and measurement equation
% x_t = C + G x_{t-1} + w_t,  w_t ~ N(0, W)
% y_t = d_t + B_t * x_t + v_t, v_t ~ N(0, V)
if LT == "GBM"

    C = [0; par.mu * deltat];
    G = [exp(-par.kappa * deltat), 0; 0, 1];
    W = [(1-exp(-2*par.kappa*deltat))/(2*par.kappa)*par.sigmachi^2, (1-exp(-par.kappa*deltat))/par.kappa * (par.sigmachi * par.sigmaxi * par.rho_chixi);...
        (1-exp(-par.kappa*deltat))/par.kappa * (par.sigmachi * par.sigmaxi * par.rho_chixi), par.sigmaxi^2 * deltat]; %the covariance matrix of w

    d1 = (1 - exp(-2 * par.kappa * ttm)) * par.sigmachi^2 / (2 * par.kappa);
    d2 = par.sigmaxi^2 * ttm;
    d3 = (1 - exp(-(par.kappa) * ttm)) * 2 * par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa);
    d = (par.mu - par.lambdaxi) * ttm - (par.lambdachi / par.kappa) * (1 - exp(-par.kappa * ttm)) + (1/2) * (d1 + d2 + d3);
    d = d';

    for i = 1:nobsn
        B1(i,:) = exp(-par.kappa * ttm(i,:));
        B2(i,:) = repelem(1, ncontracts);
        B(:,:,i) = [B1(i,:); B2(i,:)]';
    end

elseif LT == "OU"

    C = [0 ; (par.mu / par.gamma) * (1 - exp(-par.gamma * deltat))];
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
%     P0 = [100, 0; 0, 100];
    P0 = [0.01, 0.01; 0.01, 0.01];

elseif LT == "OU"
    a0 = [0; mean(y(1,:), 2)];
%     P0 = [100, 0; 0, 100];
    P0 = [1, 0; 0, 1];
%     a0 = [0 ; par.mu / par.gamma]; % a_1|0
%     P0 = [par.sigmachi^2 / (2*par.kappa), par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa + par.gamma);
%         par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa + par.gamma), par.sigmaxi^2 / (2 * par.gamma)]; % P_1_0

end

global save_ytt save_att save_ett save_vt save_att_1 save_Ptt_1 save_Ptt save_ett_1
save_ytt = zeros(nobsn, ncontracts); % Estimation of y given x_1:t, y_{t-1}
save_att = zeros(nobsn, size(a0,1)); % Updated estimation of x at time t
save_ett = zeros(nobsn, ncontracts); 
save_vt = zeros(nobsn, ncontracts); % Measurement errors
save_att_1 = zeros(nobsn, size(a0,1));
save_Ptt_1 = zeros(size(a0,1), size(a0,1), nobsn);
save_Ptt = zeros(size(a0,1), size(a0,1), nobsn);
save_ett_1 = zeros(nobsn, ncontracts); % Prediction error

att_1 = a0;
Ptt_1 = P0;
% att_1 = C + G*a0;
% Ptt_1 = G*P0*G' + W;
eLe = 0;
dLtt_1 = 0;
ett = zeros(ncontracts, 1);
% vt = 0;

for i = 1:nobsn


    %Prediction error and covariance matrix
    % e_t = y_t - E[y_t | I_{t-1}]
    ytt_1 = d(:,i) + B(:,:,i) * att_1; % E[y_t | I_{t-1}]
    yt = y(i,:)'; % y_t
    et = yt - ytt_1;

    % L_t|{t-1} = B_t P_t|{t-1} B_t' + V
    Ltt_1 = B(:,:,i) * Ptt_1 * B(:,:,i)' + V;
    %  if sum(sum(diag(eig((Ltt_1 + Ltt_1') / 2)))) > 0
    %      disp('matrix is not postive semi-definite');
    %  end

    dLtt_1 = dLtt_1 + log(det(Ltt_1)); % ln (det(L_t|{t-1}))
    eLe = eLe + et' * inv(Ltt_1) * et; % e_t' * L_t|{t-1} * e_t

    % Update equation
    % Kalman gain: K_t = P_t|{t-1} B_t' (L_t|{t-1})^(-1)
    % Expectation: a_t = a_{t|t-1} + K_t e_t
    % Covariance:
    % P_t = (I - K_t B_t) * P_{t|t-1} * (I - K_t B_t)'+ K_t * V * K_t'
    Kt = Ptt_1 * B(:,:,i)' * inv(Ltt_1);
    att = att_1 + Kt * et;
    Rt = eye(size(Ptt_1,1)) - Kt * B(:,:,i);
    Ptt = Rt * Ptt_1 * Rt' + Kt * V * Kt';

    % Forecast distributions of state variables
    % a_{t+1|t} = C + G a_t
    % P_{t+1|t} = G * P_t * G' + W
    save_att_1(i,:) = att_1';
    save_Ptt_1(:,:,i) = Ptt_1;
    att_1 = C + G * att;
    Ptt_1 = G * Ptt * G' + W;

    % Estimate y, given (x_1:t, y_{t-1})
    ytt = d(:,i) + B(:,:,i) * att;
    ett = yt - ytt; % Measurement error.

    save_ytt(i,:) = ytt';
    save_att(i,:) = att';
    save_ett(i,:) = ett';
    save_Ptt(:,:,i) = Ptt;
    save_ett_1(i,:) = et';
    save_ytt_1(i,:) = ytt_1';
    % save_et(i,:) = et';
end

logL = -(1/2) * nobsn * ncontracts * log(2*pi) - (1 / 2) * dLtt_1 - (1 / 2) * eLe;
log_L = -logL;
%------------------- End of Likelihood Function ---------------------%

