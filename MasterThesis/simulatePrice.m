%% Simulation
% This code simulated futures price dataset using:
% par : [kappa, sigmachi, lambdachi, gamma, mu, sigmaxi, lambdaxi, rho, s,
% v]
% T : maturity time, yeardays: number of days in a year, and
% monthdays : number of days in a month, for calculating time to maturities
% ncontracts : number of contracts, and nobsn : number of observations
% deltat : time difference between each observation
% model_selection : Model number
% vary_ME : Estimation method of volatilities and correlations of
% measurement errors
% ttm_by_user : if the user has own time to maturities. Otherwise 0.


function [simY, x, ttm] = simulatePrice(par, T, yeardays, monthdays, ncontracts, nobsn, deltat, model_selection, vary_ME, ttm_by_user)

kappa = par(1);
sigmachi = par(2);
lambdachi = par(3);
gamma = par(4);
mu = par(5);
sigmaxi = par(6);
lambdaxi = par(7);
rho = par(8);

if any(vary_ME == [2,4])
    vsig = par(9) + par(10)*exp(par(11)*T);
elseif any(vary_ME == [1,3])
    vsig = par(9:8+ncontracts);
end

if any(vary_ME == [1,2])
    V = diag(vsig.^2);
    
elseif vary_ME == 3
    correlv = par(9+ncontracts:end);
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
    correlv = par(12:11+ncontracts);
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

for j = 1:ncontracts
    for i = 1:nobsn
        if j == 1
            ttm(i,j) = T(j) - (i)/yeardays;
            if ttm(i,1) < 0
                ttm(i,j) = ttm((i - monthdays * 12 * T(j)), j);
            end
        else
            ttm(i,j) = T(j) - (i)/yeardays;
            if i > monthdays
                ttm(i,j) = ttm(i - monthdays, j);
            end
        end
    end
end

if ttm_by_user ~= 0 % If the time to maturity is given by the user,
    ttm = ttm_by_user;
end

if model_selection == 1
    kappa = 0;
    sigmachi = 0;
    lambdachi = 0;
    rho = 0;
    
    C = [0; (mu - (1/2) * sigmaxi^2)*deltat];
    G = [0, 0; 0, 1];
    W = [0, 0; 0, sigmaxi^2 * deltat];
    
    a1 = 0;
    a2 = sigmaxi^2 * ttm;
    a3 = 0;
    d = (mu - lambdaxi) * ttm + (1/2) * (a1 + a2 + a3);
    d = d';
    
    x(1,:) = [0; mu];
    P0 = [0, 0; 0, 0.01];
    
elseif model_selection == 2
    
    C = [0 ; (mu) * deltat];
    G = [exp(-kappa * deltat), 0; 0, 1];
    W = [(1-exp(-2*kappa*deltat))/(2*kappa)*sigmachi^2, (1-exp(-kappa*deltat))/kappa * (sigmachi * sigmaxi * rho);...
        (1-exp(-kappa*deltat))/kappa * (sigmachi * sigmaxi * rho), sigmaxi^2 * deltat]; %the covariance matrix of w
    
    x(1,:) = [0; mu];
    P0 = [0.01, 0.01; 0.01, 0.01];
    
    a1 = (1 - exp(-2 * kappa * ttm)) * sigmachi^2 / (2 * kappa);
    a2 = sigmaxi^2 * ttm;
    a3 = (1 - exp(-kappa * ttm)) * 2 * sigmachi * sigmaxi * rho / kappa;
    d = (mu - lambdaxi) * ttm - (lambdachi/kappa) * (1 - exp(-kappa * ttm)) + (1/2) * (a1 + a2 + a3);
    d = d';
    
    for i = 1:nobsn
        F1(i,:) = exp(-kappa * ttm(i,:));
    end
    
elseif model_selection == 3
    
    C = [0 ; (mu/gamma)*(1-exp(-gamma * deltat))];
    G = [exp(-kappa * deltat), 0; 0, exp(-gamma * deltat)];
    W = [(1-exp(-2*kappa*deltat))/(2*kappa)*sigmachi^2, (1-exp(-(kappa + gamma)*deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho);...
        (1-exp(-(kappa + gamma)*deltat)) / (kappa + gamma) * (sigmachi * sigmaxi * rho), (1-exp(-2*gamma*deltat))/(2*gamma)*sigmaxi^2];
    
    x(1,:) = [0 ; mu]';
    
    P0 = [sigmachi^2/2/kappa, sigmachi*sigmaxi.*rho / (kappa + gamma);
        sigmachi*sigmaxi*rho/(kappa + gamma), sigmaxi.^2/2/gamma];
    
    a1 = (1 - exp(-2 * kappa * ttm)) * sigmachi^2 / (2 * kappa);
    a2 = (1 - exp(-2 * gamma * ttm)) * sigmaxi^2 / (2 * gamma);
    a3 = (1 - exp(- (kappa + gamma) * ttm)) * 2 * sigmachi * sigmaxi * rho / (kappa + gamma);
    d = (mu - lambdaxi) / gamma * (1 - exp(-gamma * ttm)) - (lambdachi/kappa) * (1 - exp(-kappa * ttm)) + (1/2) * (a1 + a2 + a3);
    d = d';
    
    for i = 1:nobsn
        F1(i,:) = exp(-kappa * ttm(i,:));
        F2(i,:) = exp(-gamma * ttm(i,:));
        F(:, :, i) = [F1(i,:); F2(i,:)];
    end
    
else
    error("Please choose the correct model_selection option")    
end

normmu = repelem(0,size(x,2));
normsigma = [1, rho; rho, 1];
burn_in = 100;

chi(1) = x(1,1); xi(1) = x(1,2);

for i = 2:nobsn+burn_in+1
    z(i-1,:) = mvnrnd(normmu, normsigma);
    if kappa == 0
        chi(i) = 0;
    else
        chi(i) = chi(i-1) * exp(-kappa * deltat) + sigmachi * sqrt((1-exp(-2*kappa*deltat))/(2*kappa))*z(i-1,1);
    end
    if gamma == 0
        xi(i) = xi(i-1) + mu * deltat  + sigmaxi * sqrt(deltat)*z(i-1,2);
    else
        xi(i) = xi(i-1) * exp(-gamma * deltat) + mu / gamma * (1-exp(-gamma * deltat)) + sigmaxi * sqrt((1-exp(-2*gamma*deltat))/(2*gamma))*z(i-1,2);
    end
end

x = [chi;xi]';
x = x(burn_in+2:end,:); %Burn-in period removed.

vnormmu = repelem(0, ncontracts);
vnormsigma = V;

for i = 1:nobsn
    v(i,:) = mvnrnd(vnormmu, vnormsigma);
    if model_selection == 1
        y(i,:) = d(:,i) + [repelem(1,ncontracts); repelem(1,ncontracts)]' * x(i,:)' + v(i,:)';
    elseif model_selection == 2
        y(i,:) = d(:,i) + [F1(i,:); repelem(1,ncontracts)]' * x(i,:)' + v(i,:)';
    elseif model_selection == 3
        y(i,:) = d(:,i) + [F1(i,:); F2(i,:)]' * x(i,:)' + v(i,:)';
    end
end

simY = y;
timetomaturity = ttm;