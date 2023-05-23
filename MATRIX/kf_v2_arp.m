%% Kalman filter

function log_L = kf_v2_arp(par_init, par_names, y, deltat, ttm, LT, correlation, serial)

% Defining the dimensions of the dataset.
[nobsn, ncontracts] = size(y);
n_par = size(par_init, 2);
par = set_parameters(LT, ncontracts, par_names, par_init, correlation, serial);
if serial == "yes"
    n_lag = size(par.phi,2)/ncontracts;
end

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
    d_temp = (par.mu - par.lambdaxi) * ttm - (par.lambdachi / par.kappa) * (1 - exp(-par.kappa * ttm)) + (1/2) * (d1 + d2 + d3);
    d_temp = d_temp';

    for i = 1:nobsn
        B1(i,:) = exp(-par.kappa * ttm(i,:));
        B2(i,:) = repelem(1, ncontracts);
        B_temp(:,:,i) = [B1(i,:); B2(i,:)]';
    end

elseif LT == "OU"

    C = [0 ; (par.mu / par.gamma) * (1 - exp(-par.gamma * deltat))];
    G = [exp(-par.kappa * deltat), 0; 0, exp(-par.gamma * deltat)];
    W = [(1 - exp(-2 * par.kappa * deltat)) / (2 * par.kappa) * par.sigmachi^2, (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi);
        (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi), (1 - exp(-2 * par.gamma * deltat)) / (2 * par.gamma) * par.sigmaxi^2];

    d1 = (1 - exp(-2 * par.kappa * ttm)) * par.sigmachi^2 / (2 * par.kappa);
    d2 = (1 - exp(-2 * par.gamma * ttm)) * par.sigmaxi^2 / (2 * par.gamma);
    d3 = (1 - exp(-(par.kappa + par.gamma) * ttm)) * 2 * par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa + par.gamma);
    d_temp = (par.mu - par.lambdaxi) / par.gamma * (1 - exp(-par.gamma * ttm)) - (par.lambdachi / par.kappa) * (1 - exp(-par.kappa * ttm)) + (1/2) * (d1 + d2 + d3);
    d_temp = d_temp';

    for i = 1:nobsn
        B1(i,:) = exp(-par.kappa * ttm(i,:));
        B2(i,:) = exp(-par.gamma * ttm(i,:));
        B_temp(:,:,i) = [B1(i,:); B2(i,:)]';
    end

else
    error('Please specify the process of the long-term factor LT.')
end

vsig2 = par.s.^2;
if correlation == 0
    V_temp = diag(vsig2);
    
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
    V_temp = D^(1/2) * CorMat * D^(1/2);

else
    error('correlation must be 0 or 1.')
end

d = []; B = []; Ct = [];

for t = 1+n_lag:nobsn
    phiy = zeros(ncontracts, 1); phiBG = zeros(ncontracts,1);
    phid = zeros(ncontracts, 1); Ginv = zeros(size(G,1), size(G,2));
    phiBGinv = zeros(ncontracts, size(G,2)); phiBGWG = zeros(ncontracts,ncontracts);
    GWG = zeros(size(Ginv,1), size(G,2)); GBphi = zeros(size(G,1), ncontracts);
    for j = 1:n_lag
        phiy = phiy + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*y(t-j,:)';
        phid = phid + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*d_temp(:,t-j);
        for k = 1:j
            Ginv = Ginv + inv(G)^k;
        end
        phiBG = phiBG + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*B_temp(:,:,t-j)*Ginv*C;
        phiBGinv = phiBGinv + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*B_temp(:,:,t-j)*inv(G)^j;
        GBphi = GBphi + (inv(G)')^(j-1)*B_temp(:,:,t-j)'*diag(par.phi((j-1)*ncontracts+1:j*ncontracts));
    end
    save_phiBGinv(:,:,t) = phiBGinv;
    d(:,t) = phiy + d_temp(:,t) - phid + phiBG;
    B(:,:,t) = B_temp(:,:,t) - phiBGinv;
    Ct(:,:,t) = W*inv(G)'*GBphi;
end

% Kalman Filter
% Initial distribution of state variables
% E(x_{1|0}) = a0
% Cov(x_{1|0}) = P0
if LT == "GBM"
    a0 = [0; y(1,1)];
    P0 = [0.01, 0.01; 0.01, 0.01];

elseif LT == "OU"
    a0 = [0; y(1,1)];
    P0 = [0.01, 0.01; 0.01, 0.01];

end

global save_ytt save_att save_ett save_vtt
save_ytt = zeros(nobsn, ncontracts); % Estimation of y given x_1:t, y_{t-1}
save_att = zeros(nobsn, length(a0)); % Updated estimation of x at time t
save_ett = zeros(nobsn, ncontracts); % Prediction error
save_vtt = zeros(nobsn, ncontracts); % Measurement errors

att_1(:,1) = a0;
Ptt_1(:,:,1) = P0;

for i = 1:n_lag
    ModX = fitlm(B_temp(:,:,i), y(i,:)' - d_temp(:,i), 'intercept', false);
    att(:,i) = ModX.Coefficients{:,i};
    K = Ptt_1(:,:,i)*B_temp(:,:,i)'*(inv(B_temp(:,:,i)*Ptt_1(:,:,i)*B_temp(:,:,i)'+V_temp));
    Ptt(:,:,i) = (eye(2)-K*B_temp(:,:,i))*Ptt_1(:,:,i);
    att_1(:,i+1) = C + G * att(:,i);
    Ptt_1(:,:,i+1)  = G * Ptt(:,:,i) * G' + W;
end

eLe = 0;
dLtt_1 = 0;
vtt = zeros(nobsn, ncontracts)';

for i = 1+n_lag:nobsn

if n_lag == 1
    vt = zeros(ncontracts, nobsn);

elseif n_lag > 1
    Gw = zeros(length(C),1);
    phiBGw = zeros(ncontracts, 1);
    for r = 2:n_lag
        for k = 2:r
            wt = ks_v1(i, k, LT, par, deltat, ttm, att, att_1, Ptt, Ptt_1);
            Gw = Gw + inv(G)^(r-k+1) * wt(:,k); %% Need to fix this section.
        end
        phiBGw = phiBGw + diag(par.phi((r-1)*ncontracts+1:r*ncontracts)) * B_temp(:,:,i-2) * Gw;
    end
    vt(:,i) = phiBGw;

end

    %Prediction error and covariance matrix
    % e_t = y_t - E[y_t | I_{t-1}]

    ytt_1 = d(:,i) + B(:,:,i) * att_1(:,i) + vt(:,i); % E[y_t | I_{t-1}]
    yt = y(i,:)'; % y_t
    et = yt - ytt_1;
    
    phiBGinv = zeros(ncontracts, size(G,2));
    for j = 1:n_lag
        phiBGinv = phiBGinv + diag(par.phi((j-1)*ncontracts+1:j*ncontracts))*B_temp(:,:,i-j)*inv(G)^j;
    end
    V = phiBGinv * W * phiBGinv' + V_temp;
    Ltt_1 = B(:,:,i) * Ptt_1(:,:,i) * B(:,:,i)' + V + B(:,:,i)*W*phiBGinv' + phiBGinv * W' * B(:,:,i)';
%      if sum(sum(diag(eig((Ltt_1 + Ltt_1') / 2)))) > 0
%          disp('matrix is not postive semi-definite');
%      end

    dLtt_1 = dLtt_1 + log(det(Ltt_1)); % ln (det(L_t|{t-1}))
    eLe = eLe + et' * inv(Ltt_1) * et; % e_t' * L_t|{t-1} * e_t

    % Update equation
    % Kalman gain: K_t = P_t|{t-1} B_t' (L_t|{t-1})^(-1)
    % Expectation: a_t = a_{t|t-1} + K_t e_t
    % Covariance:
    % P_t = (I - K_t B_t) * P_{t|t-1} * (I - K_t B_t)'+ K_t * V * K_t'
    Kt = (Ptt_1(:,:,i) * B(:,:,i)' + Ct(:,:,i)) * inv(Ltt_1 + B(:,:,i-1) * Ct(:,:,i) + Ct(:,:,i)'*B(:,:,i)');
    att(:,i) = att_1(:,i) + Kt * et;
    Rt = eye(size(Ptt_1(:,:,i),1)) - Kt * B(:,:,i);
    Ptt(:,:,i) = Rt * Ptt_1(:,:,i) * Rt' - Rt * Ct(:,:,i) * Kt' - Kt * Ct(:,:,i)' * Rt' + Kt * V * Kt';

    % Forecast distributions of state variables
    % a_{t+1|t} = C + G a_t
    % P_{t+1|t} = G * P_t * G' + W
    att_1(:,i+1) = C + G * att(:,i);
    Ptt_1(:,:,i+1) = G * Ptt(:,:,i) * G' + W;
    
    J = Ptt(:,:,i-1) * G' * inv(Ptt_1(:,:,i));
    xt_1t = att(:,i-1) + J * (att(:,i)-att_1(:,i));
    wtt = att(:,i) - C - G * xt_1t;
    ytt = d(:,i) + B(:,:,i) * att(:,i) + diag(par.phi)*B_temp(:,:,i-1)*inv(G)*wtt;
    ett = yt - ytt; % Measurement error.
    vtt(:,i) = diag(par.phi)*vtt(:,i-1) + ett;
    
    save_ytt(i,:) = ytt';
    save_att(i,:) = att(:,i)';
    save_ett(i,:) = ett';
    save_dLtt_1(:,i) = dLtt_1;
end

save_ytt(1:n_lag,:) = y(1:n_lag,:);
save_att(1:n_lag, :) = att(:,1:n_lag)';
save_vtt = vtt';
logL = - (1 / 2) * dLtt_1 - (1 / 2) * eLe;
log_L = -logL;
%------------------- End of Likelihood Function ---------------------%

