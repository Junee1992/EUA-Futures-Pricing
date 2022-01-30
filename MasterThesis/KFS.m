%% Kalman Filter and Smoothing

% This code produces the estimated values of latent variables (bivariate) using Kalman
% Filter and the Kalman Smoothing in the linear state-space model, and also
% calculates the log likelihood.

function log_L = KFS(par, y, deltat, ttm_by_user, model_selection, vary_ME, mat)

ncontracts = size(y,2);
nobsn = size(y,1);

kappa = par(1);
sigmachi = par(2);
lambdachi = par(3);
gamma = par(4);
mu = par(5);
sigmaxi = par(6);
lambdaxi = par(7);
rho = par(8);

if any(model_selection == [1,2])
    ttm = ttm_by_user;
    a0 = [0; (y(1,1))]; % Uses the first y-value as the initial value for the long-term factor.
    
    % The state equation
    % x(t) = C + G * x(t-1) + w(t), w(t) ~ N(0,W)
    
    if model_selection == 1
        P0 = [0, 0; 0, 0.01];
        C = [0; (mu - (1/2) * sigmaxi^2)*deltat];
        G = [0, 0; 0, 1];
        W = [0, 0; 0, sigmaxi^2 * deltat];
        
        a1 = 0;
        a2 = sigmaxi^2 * ttm;
        a3 = 0;
        d = (mu - lambdaxi) * ttm + (1/2) * (a1 + a2 + a3);
        d = d';
        
        F1 = zeros(nobsn, ncontracts);
        
    elseif model_selection == 2
        P0 = [0.01, 0.01; 0.01, 0.01];
        C = [-(lambdachi/kappa) * (1 - exp(-kappa * deltat)); (mu - lambdaxi) * deltat];
        G = [exp(-kappa * deltat), 0; 0, 1];
        W = [(1-exp(-2*kappa*deltat))/(2*kappa)*sigmachi^2, (1-exp(-kappa*deltat))/kappa * (sigmachi * sigmaxi * rho);...
            (1-exp(-kappa*deltat))/kappa * (sigmachi * sigmaxi * rho), sigmaxi^2 * deltat]; %the covariance matrix of w
        
        a1 = (1 - exp(-2 * kappa * ttm)) * sigmachi^2 / (2 * kappa);
        a2 = sigmaxi^2 * ttm;
        a3 = (1 - exp(-kappa * ttm)) * 2 * sigmachi * sigmaxi * rho / kappa;
        d = (mu - lambdaxi) * ttm - (lambdachi/kappa) * (1 - exp(-kappa * ttm)) + (1/2) * (a1 + a2 + a3);
        d = d';
        
        for i = 1:nobsn
            F1(i,:) = exp(-kappa * ttm(i,:));
        end
    end
    
elseif model_selection == 3
    ttm = ttm_by_user;
    
    a0 = [0 ; mu / gamma];
    P0 = [sigmachi^2/2/kappa, sigmachi*sigmaxi*rho / (kappa + gamma);
        sigmachi*sigmaxi*rho/(kappa + gamma), sigmaxi^2/2/gamma]; % P_1|0
    
    C = [0; (mu / gamma) * (1 - exp(-gamma * deltat))];
    G = [exp(-kappa * deltat), 0; 0, exp(-gamma * deltat)];
    W = [(1-exp(-2*kappa*deltat))/(2*kappa)*sigmachi^2, (1-exp(-(kappa + gamma)*deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho);...
        (1-exp(-(kappa + gamma)*deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho), (1-exp(-2*gamma*deltat))/(2*gamma)*sigmaxi^2]; %the covariance matrix of w
    
    a1 = (1 - exp(-2 * kappa * ttm)) * sigmachi^2 / (2 * kappa);
    a2 = (1 - exp(-2 * gamma * ttm)) * sigmaxi^2 / (2 * gamma);
    a3 = (1 - exp(- (kappa + gamma) * ttm)) * 2 * sigmachi * sigmaxi * rho / (kappa + gamma);
    d = (mu - lambdaxi) / gamma * (1 - exp(-gamma * ttm)) - (lambdachi/kappa) * (1 - exp(-kappa * ttm)) + (1/2) * (a1 + a2 + a3);
    d = d';
    
    for i = 1:nobsn
        F1(i,:) = exp(-kappa * ttm(i,:));
        F2(i,:) = exp(-gamma * ttm(i,:));
    end

else
    error("Please choose the correct vary_ME option")    
end

if vary_ME == 1
    vsig = par(9:8+ncontracts);
    V = diag(vsig.^2);
    
elseif vary_ME == 2
    f = par(9);
    g = par(10);
    h = par(11);
    
    for i = 1:ncontracts
        vsig(i) = f + g * exp(-h * mat(i));
    end
    V = diag(vsig.^2);
    
elseif vary_ME == 3
    vsig = par(9:8+ncontracts);
    correlv = par(9+ncontracts:8+2*ncontracts);
    CorMat = diag(repelem(1,ncontracts));
    
    for i = 1:ncontracts
        for j = 1:ncontracts
            if i == j
                CorMat(i,j) = 1;
            else
                CorMat(i,j) = correlv(i) * correlv(j);
            end
        end
    end
    
    D = diag(vsig.^2);
    V = D.^(1/2) * CorMat * D.^(1/2);
    
elseif vary_ME == 4
    f = par(9);
    g = par(10);
    h = par(11);
    
    for i = 1:ncontracts
        vsig(i) = f + g * exp(-h * mat(i));
    end
    
    correlv = par(12:end);
    CorMat = diag(repelem(1,ncontracts));
    
    for i = 1:ncontracts
        for j = 1:ncontracts
            if i == j
                CorMat(i,j) = 1;
            else
                CorMat(i,j) = correlv(i) * correlv(j);
            end
        end
    end
    
    D = diag(vsig.^2);
    V = D.^(1/2) * CorMat * D.^(1/2);
    
end

% Kalman Filter
att_1 = a0;
Ptt_1 = P0;

global save_xf save_ett save_ytt save_xs_MBF save_xs_RTS

save_ett = zeros(nobsn, ncontracts);
save_et = zeros(nobsn, ncontracts);
save_Ptt_1 = [];
save_Ptt = [];
save_Ltt_1 = zeros(ncontracts, ncontracts, nobsn);
save_dLtt_1 = zeros(nobsn, 1);
save_eLe = zeros(nobsn,1);
save_ytt = zeros(nobsn, ncontracts);
save_att = zeros(nobsn,size(a0,1));
save_att_1 = zeros(nobsn, size(a0,1));
save_Ft = zeros(size(a0,1), ncontracts, nobsn);
save_att_1(1, :) = a0';
save_Ptt_1(:, :, 1) = P0;
save_Kt = zeros(size(a0,1), ncontracts, nobsn);

if any(model_selection == [1,2])
    
    for i = 1:nobsn
        
        %Prediction error and covariance matrix
        Ft = [F1(i,:); repelem(1, ncontracts)]';
        save_Ft(:, :, i) = Ft';
        ytt_1 = d(:,i) + [F1(i,:); repelem(1, ncontracts)]'*att_1; %yt|t-1
        yt = y(i,:)';
        et = yt - ytt_1; %et|t-1
        
        Ltt_1 = [F1(i,:); repelem(1, ncontracts)]' * Ptt_1 * [F1(i,:); repelem(1, ncontracts)] + V; %Lt|t-1
        dLtt_1 = det(Ltt_1); %det(Lt|t-1)
        
        %Updating equation
        Kt = Ptt_1 * [F1(i,:); repelem(1, ncontracts)] * inv(Ltt_1);
        att = att_1 + Kt * et; %at|t
        Ptt = Ptt_1 - Kt * [F1(i,:); repelem(1, ncontracts)]' * Ptt_1; %Pt|t
        
        att_1 = C + G * att; %at|t-1
        Ptt_1 = G * Ptt * G' + W; %Pt|t-1
        
        % The measurement equation
        ytt = [F1(i,:); repelem(1, ncontracts)]' * att + d(:,i); %yhat
        ett = yt - ytt; %et|t
        
        save_ytt(i,:) = ytt';
        save_ett(i,:) = ett';
        save_et(i,:) = et';
        save_att(i,:) = att';
        save_att_1(i+1,:) = att_1';
        save_Ptt_1(:, :, i+1) = Ptt_1;
        save_Ptt(:, :, i) = Ptt;
        save_dLtt_1(i,:) = dLtt_1;
        save_Ltt_1(:, :, i) = Ltt_1;
        save_eLe(i,:) = et' * inv(Ltt_1) * et;
        save_Ft(:, :, i) = Ft';
        save_Kt(:, :, i) = Kt;
        
    end
    
elseif model_selection == 3
    
    for i = 1:nobsn
        
        %Prediction error and covariance matrix
        Ft = [F1(i,:); F2(i,:)]';
        ytt_1 = d(:,i) + [F1(i,:); F2(i,:)]'*att_1;
        yt = y(i,:)';
        et = yt - ytt_1; %et|t-1
        
        Ltt_1 = [F1(i,:); F2(i,:)]' * Ptt_1 * [F1(i,:); F2(i,:)] + V; %Lt|t-1
        
        if  sum(sum(diag(eig((Ltt_1+Ltt_1')/2)<0)))>0
            disp('matrix is not positive semi-definite');
        end
        
        dLtt_1 = det(Ltt_1); %det(Lt|t-1)
        
        %Updating equation
        Kt = Ptt_1 * [F1(i,:); F2(i,:)] * inv(Ltt_1);
        att = att_1 + Kt * et; %at|t
        Ptt = Ptt_1 - Kt * [F1(i,:); F2(i,:)]' * Ptt_1; %Pt|t
        
        att_1 = C + G * att; %at|t-1
        Ptt_1 = G * Ptt * G' + W; %Pt|t-1
        
        ytt = [F1(i,:); F2(i,:)]' * att + d(:,i);
        ett = yt - ytt; %et|t
        
        save_ytt(i,:) = ytt';
        save_ett(i,:) = ett';
        save_et(i,:) = et';
        save_att(i,:) = att';
        save_att_1(i+1,:) = att_1';
        save_Ptt_1(:, :, i+1) = Ptt_1;
        save_Ptt(:, :, i) = Ptt;
        save_dLtt_1(i,:) = dLtt_1;
        save_Ltt_1(:, :, i) = Ltt_1;
        save_eLe(i,:) = et' * inv(Ltt_1) * et;
        save_Ft(:, :, i) = Ft';
        save_Kt(:, :, i) = Kt;
        
    end
end

save_xf = save_att;

% Kalman Smoother

% MBF Smoother
save_xs_MBF = zeros(nobsn, size(a0,1));
save_xs_MBF(nobsn,:) = save_att(nobsn,:);
save_Ps_MBF = zeros(size(a0,1), size(a0,1), nobsn);
save_Ps_MBF(:, :, nobsn) = save_Ptt(:, :, nobsn);

r = zeros(size(a0,1), 1, nobsn);
R = zeros(size(a0,1), size(a0,1), nobsn);

for i = nobsn : -1 : 1
    if i == nobsn
        r = [0; 0];
        R = [0, 0; 0, 0];
    end
    
    r = save_Ft(:, :, i) * inv(save_Ltt_1(:, :, i)) * save_et(i,:)' + (G - G * save_Kt(:, :, i) * save_Ft(:, :, i)')' * r;
    R = save_Ft(:, :, i) * inv(save_Ltt_1(:, :, i)) * save_Ft(:, :, i)' + (G - G * save_Kt(:, :, i) * save_Ft(:, :, i)')' * ...
        R * (G - G * save_Kt(:, :, i) * save_Ft(:, :, i)');
    
    xs = save_att_1(i,:)' + save_Ptt_1(:, :, i) * r;
    Ps = save_Ptt_1(:, :, i) - save_Ptt_1(:, :, i) * R * save_Ptt_1(:, :, i);
    
    save_xs_MBF(i,:) = xs';
    save_Ps_MBF(:, :, i) = Ps;
end
    

% RTS smoother
save_xs_RTS = zeros(nobsn, size(a0,1));
save_xs_RTS(nobsn,:) = save_att(nobsn,:);
save_Ps_RTS = zeros(size(a0,1), size(a0,1), nobsn);
save_Ps_RTS(:, :, nobsn) = save_Ptt(:, :, nobsn);

if model_selection == 1
    save_xs_RTS = 0;
else
    for i = nobsn : -1 : 2
        J = save_Ptt(:, :, i-1) * G' * inv(save_Ptt_1(:, :, i));
        xs = save_att(i-1, :)' + J * (save_xs_RTS(i, :)' - save_att_1(i, :)');
        Ps = save_Ptt(:, :, i-1) + J * (save_Ps_RTS(:, :, i) - save_Ptt_1(:, :, i)) * J';
        save_xs_RTS(i-1,:) = xs';
        save_Ps_RTS(:, :, i-1) = Ps;
    end
end

logL = -(nobsn * ncontracts/2) * log(2*pi) - (1/2)*sum(log(save_dLtt_1)) - (1/2) * sum(save_eLe);
log_L = -logL;

    
