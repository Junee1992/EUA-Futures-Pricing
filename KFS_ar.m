%% Kalman Filter and Parameter Estimation of Schwartz-Smith Two-Factor Model
% This code produces the estimated values of latent variables using
% Kalman filter in the linear state-space model, and also computes
% the negative log-likelihood function for parameter estimation.

% This code can incorporate models with serially correlated measurement errors,
% where those measurement errors follow an AR(1) process.

% In mathematical expression: y_t = d_t + B_t x_t + v_t,
% where v_t = M_1 v_{t-1} + epsilon_t, epsilon_t ~ N(0,V) %

% This code only allows for two OU processes for two state variables.

function [log_L] = KFS_ar(par, y, deltat, ttm, model, a0, P0, ms_error, AR_process)

% Defining the dimension of our dataset, and time to maturity.
ncontracts = size(y,2);
nobsn = size(y,1);

if size(par, 2) < 8 + ncontracts
    error('Need at least (8+ncontracts) parameters.')
end

% Parameters
kappa = par(1);
sigmachi = par(2);
lambdachi = par(3);
gamma = par(4);
mu = par(5);
sigmaxi = par(6);
lambdaxi = par(7);
rho = par(8);

if model == "price"
    
    % Initial distribution of state variables
    % E(x_{1|0}) = a0
    % Cov(x_{1|0}) = P0
    a0 = [0 ; mu / gamma]; % a_1|0
    P0 = [sigmachi^2 / (2*kappa), sigmachi * sigmaxi * rho / (kappa + gamma);
        sigmachi * sigmaxi * rho / (kappa + gamma), sigmaxi^2 / (2 * gamma)]; % P_1_0
    
    % Transition matrices in state equation
    % x_t = C + G x_{t-1} + w_t, w_t ~ N(0, W)
    C = [0 ; (mu / gamma) * (1 - exp(-gamma * deltat))];
    G = [exp(-kappa * deltat), 0;
        0, exp(-gamma * deltat)];
    W = [(1 - exp(-2 * kappa * deltat)) / (2 * kappa) * sigmachi^2, (1 - exp(-(kappa + gamma) * deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho);
        (1 - exp(-(kappa + gamma) * deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho), (1 - exp(-2 * gamma * deltat)) / (2 * gamma) * sigmaxi^2];
    
    % Constant and transition matrices in measurement equation
    % y_t = d_t + B_t x_t + v_t + M_1 v_{t-1}, v ~ N(0, V)
    d1 = (1 - exp(-2 * kappa * ttm)) * sigmachi^2 / (2 * kappa);
    d2 = (1 - exp(-2 * gamma * ttm)) * sigmaxi^2 / (2 * gamma);
    d3 = (1 - exp(-(kappa + gamma) * ttm)) * 2 * sigmachi * sigmaxi * rho / (kappa + gamma);
    d = (mu - lambdaxi) / gamma * (1 - exp(-gamma * ttm)) - (lambdachi / kappa) * (1 - exp(-kappa * ttm)) + (1/2) * (d1 + d2 + d3);
    d = d';
    
    B = zeros(ncontracts, size(a0,1), nobsn);
    for i = 1:nobsn
        B1(i,:) = exp(-kappa * ttm(i,:));
        B2(i,:) = exp(-gamma * ttm(i,:));
        B(:,:,i) = [B1(i,:); B2(i,:)]';
    end
    
    % Covariance matrix of measurement error
    vsig2 = par(9:8+ncontracts).^2;
    
    if ms_error == 1
        V = diag(vsig2); % Assumes independence in volatilities of measurement errors
        
        if size(par, 2) == 8 + ncontracts
            par = par;
        else
            if AR_process == "noAR"
                error('Cannot use ms_error = 1 with "no AR". Try ms_error = 2.')
            elseif AR_process == "AR"
                AR_coeff = par(9+ncontracts:end);
                disp('Check if you have correctly entered AR coefficients, not intercorrelation coefficients.')
                disp('Stop the process if you have intercorrelation coefficients in your parameter set.')
            end
        end
        
        m_diag = diag(AR_coeff);

    elseif ms_error == 2
        correl = par(9+ncontracts:8+2*ncontracts); % Inter-correlated measurement errors
        
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
        
        if size(par, 2) == 8 + 2 * ncontracts
            par = par;
        else
            AR_coeff = par(9+2*ncontracts:end);
        end
        
        m_diag = diag(AR_coeff);
        
    end
    
    % Kalman Filter
    % Initial prior distribution of the state vector
    att_1 = a0;
    Ptt_1 = P0;
    
    global save_ytt save_att save_ett
    save_ytt = zeros(nobsn, ncontracts); % Estimation of y given x_1:t, y_{t-1}
    save_att = zeros(nobsn, size(a0,1)); % Updated estimation of x at time t
    save_ett = zeros(nobsn, ncontracts); % Prediction error
    
    eLe = 0;
    dLtt_1 = 0;
    
    for i = 1:nobsn
        
        %Prediction error and covariance matrix
        % e_t = y_t - E[y_t | I_{t-1}]
        ytt_1 = d(:,i) + B(:,:,i) * att_1; % E[y_t | I_{t-1}]
        yt = y(i,:)'; % y_t
        
        if exist('m_diag')
            if i == 1
                et = yt - ytt_1;
            else
                et = yt - ytt_1 - m_diag * ett;
            end
        else
            et = yt - ytt_1;
        end
        
        % L_t|{t-1} = B_t P_t|{t-1} B_t' + V
        Ltt_1 = B(:,:,i) * Ptt_1 * B(:,:,i)' + V;
        
        dLtt_1 = dLtt_1 + log(det(Ltt_1)); % ln (det(L_t|{t-1}))
        eLe = eLe + et' * inv(Ltt_1) * et; % e_t' * L_t|{t-1} * e_t
        
        % Update equation
        % Kalman gain: K_t = P_t|{t-1} B_t' (L_t|{t-1})^(-1)
        % Expectation: a_t = a_{t|t-1} + K_t e_t
        % Covariance: P_t = (I - K_t B_t)' * P_{t|t-1}
        Kt = Ptt_1 * B(:,:,i)' * inv(Ltt_1);
        att = att_1 + Kt * et;
        Ptt = Ptt_1 - Kt * B(:,:,i) * Ptt_1;
        
        % Forecast distributions of state variables
        % a_{t+1|t} = C + G a_t
        % P_{t+1|t} = G * P_t * G' + W
        att_1 = C + G * att;
        Ptt_1 = G * Ptt * G' + W;
        
        % Estimate y, given (x_1:t, y_{t-1})
        ytt = d(:,i) + B(:,:,i) * att;
        ett = yt - ytt; % Measurement error.
        
        save_ytt(i,:) = ytt';
        save_att(i,:) = att';
        save_ett(i,:) = ett';
        
    end
    
    logL = - (1 / 2) * dLtt_1 - (1 / 2) * eLe;
    log_L = -logL;
    
elseif model == "return"
    [nobsn, ncontracts] = size(y);
    
    kappa = par(1);
    sigmachi = par(2);
    lambdachi = par(3);
    gamma = par(4);
    mu = par(5);
    sigmaxi = par(6);
    lambdaxi = par(7);
    rho = par(8);
    
    C = [0; (mu / gamma) * (1 - exp(-gamma * deltat)); 0 ; 0];
    G = [exp(-kappa * deltat), 0, 0, 0; 0, exp(-gamma * deltat), 0, 0; 1, 0, 0, 0; 0, 1, 0, 0];
    W = [(1-exp(-2*kappa*deltat))/(2*kappa)*sigmachi^2, (1-exp(-(kappa + gamma)*deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho), 0, 0;...
        (1-exp(-(kappa + gamma)*deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho), (1-exp(-2*gamma*deltat))/(2*gamma)*sigmaxi^2, 0, 0;...
        0, 0, 0, 0; 0, 0, 0, 0]; %the covariance matrix of w
    
    a1 = (1 - exp(-2 * kappa * ttm)) * sigmachi^2 / (2 * kappa);
    a2 = (1 - exp(-2 * gamma * ttm)) * sigmaxi^2 / (2 * gamma);
    a3 = (1 - exp(- (kappa + gamma) * ttm)) * 2 * sigmachi * sigmaxi * rho / (kappa + gamma);
    d_old = (mu - lambdaxi) / gamma * (1 - exp(-gamma * ttm)) - (lambdachi/kappa) * (1 - exp(-kappa * ttm)) + (1/2) * (a1 + a2 + a3);
    d_old = d_old';
    for i = 1:nobsn-1
        d(:,i) = d_old(:,i+1) - d_old(:,i);
    end
    
    for i = 2:nobsn
        F1(i-1,:) = exp(-kappa * ttm(i,:));
        F2(i-1,:) = exp(-gamma * ttm(i,:));
        F3(i-1,:) = -exp(-kappa * ttm(i-1,:));
        F4(i-1,:) = -exp(-gamma * ttm(i-1,:));
        F(:,:,i-1) = [F1(i-1,:); F2(i-1,:); F3(i-1,:); F4(i-1,:)];
    end
    
    % Covariance matrix of measurement error
    vsig2 = par(9:8+ncontracts).^2;
    
    if ms_error == 1
        V = diag(vsig2); % Assumes independence in volatilities of measurement errors
        
        if size(par, 2) == 8 + ncontracts
            par = par;
        else
            if AR_process == "noAR"
                error('Cannot use ms_error = 1 with "no AR". Try ms_error = 2.')
            elseif AR_process == "AR"
                AR_coeff = par(9+ncontracts:end);
                disp('Check if you have correctly entered AR coefficients, not intercorrelation coefficients.')
                disp('Stop the process if you have intercorrelation coefficients in your parameter set.')
            end
        end
        
        m_diag = diag(AR_coeff);
         
    elseif ms_error == 2
        correl = par(9+ncontracts:8+2*ncontracts); % Inter-correlated measurement errors
        
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
        
        if size(par, 2) == 8 + 2 * ncontracts
            par = par;
        else
            AR_coeff = par(9+2*ncontracts:end);
        end
        
        m_diag = diag(AR_coeff);
        
    end
    
    % Kalman Filter
    att_1 = a0;
    Ptt_1 = P0;
    
    global save_ytt save_att save_ett
    save_ytt = zeros(nobsn, ncontracts); % Estimation of y given x_1:t, y_{t-1}
    save_att = zeros(nobsn, size(a0,1)); % Updated estimation of x at time t
    save_ett = zeros(nobsn, ncontracts); % Prediction error
     
    y_return = diff(y);
    
    for i = 1:nobsn-1
        
        %Prediction error and covariance matrix
        ytt_1 = d(:,i) + [F(:,:,i)]'*att_1;
        yt = y_return(i,:)';
        
        if exist('m_diag')
            if i == 1
                et = yt - ytt_1;
            else
                et = yt - ytt_1 - m_diag * ett;
            end
        else
            et = yt - ytt_1;
        end
        
        Ltt_1 = [F(:,:,i)]' * Ptt_1 * [F(:,:,i)] + V; %Lt|t-1
         
        dLtt_1 = det(Ltt_1); %det(Lt|t-1)
        
        %Updating equation
        Kt = Ptt_1 * [F(:,:,i)] * inv(Ltt_1);
        att = att_1 + Kt * et; %at|t
        Ptt = Ptt_1 - Kt * [F(:,:,i)]' * Ptt_1; %Pt|t
        
        att_1 = C + G * att; %at|t-1
        Ptt_1 = G * Ptt * G' + W; %Pt|t-1
        
        ytt = [F(:,:,i)]' * att + d(:,i);
        ett = yt - ytt; %et|t
        
        save_ytt(i,:) = ytt';
        save_ett(i,:) = ett';
        save_et(i,:) = et';
        save_att(i,:) = att';
        save_dLtt_1(i,:) = dLtt_1;
        save_eLe(i,:) = et' * inv(Ltt_1) * et;
    end
    
    logL = - (1/2)*sum(log(save_dLtt_1)) - (1/2) * sum(save_eLe);
    log_L = -logL;
end

%--------------------- End of Likelihood Function -----------------------%
