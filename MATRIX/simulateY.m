%% Simulation for futures prices with serially correlated measurement errors
% The simulation is done under the setting of Schwartz-Smith Two-Factor
% model, with two state variables following O-U processes.
% Also, we assume the measurement errors are serially correlated with AR(1)
% process.

% Mathematical expressions for the linear state-space model
% x_t = C + G * x_{t-1}+ + w_t, w_t ~ N(0, W)
% y_t = d_t + B_t * x_{t-1} + M * v_{t-1} + e_t, v_t ~ N(0, V) and
% e_t ~ N(0,V_e)

% For simplicity, we assume 260 days in a year, 30 days in a month.

function [simY, x, ttm, v] = simulateY(par, par_names, T, ncontracts, nobsn, deltat, LT, correlation, serial)
% par is a vector of relevant parameters of the model.
% T is a vector of maturity time in years. For example, 1/12 = 1 month.
% ncontracts is the number of contracts, and nobsn is the number of
% observations.
% deltat is the time difference between each observation.
%i.e. deltat = t_i - t_{i-1}.

n_par = size(par, 2);
par = set_parameters(LT, ncontracts, par_names, par, correlation, serial);
if serial == "yes"
    n_lags = length(par.phi)/ncontracts;
else
    n_lags = 0;
end

% Time to maturity
ttm = zeros(nobsn, ncontracts);
yeardays = 260;
monthdays = 30;

for j = 1:ncontracts
    for i = 1:nobsn
        if j == 1
            ttm(i,j) = T(j) - i/yeardays;
            if ttm(i,j) < 0
%                 ttm(i,j) = ttm((i - monthdays * 12 * T(j)), j);
%                 ttm(i,j) = ttm(i-yeardays/2,j);
                ttm(i,j) = ttm(i-yeardays,j);
            end
        else
            ttm(i,j) = T(j) - i/yeardays;
%             if i > monthdays
%             if i > yeardays/2
                if i > yeardays
%                 ttm(i,j) = ttm(i - yeardays/2, j);
                ttm(i,j) = ttm(i - yeardays, j);
            end
        end
    end
end


if LT == "GBM"

    C = [0; par.mu * deltat];
    G = [exp(-par.kappa * deltat), 0; 0, 1];
    W = [(1-exp(-2*par.kappa*deltat))/(2*par.kappa)*par.sigmachi^2, (1-exp(-par.kappa*deltat))/par.kappa * (par.sigmachi * par.sigmaxi * par.rho_chixi);...
        (1-exp(-par.kappa*deltat))/par.kappa * (par.sigmachi * par.sigmaxi * par.rho_chixi), par.sigmaxi^2 * deltat]; %the covariance matrix of w

    x = zeros(nobsn, size(C,1));
    x(1,:) = [0 ; par.mu]'; % Initial value of x.

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

    x = zeros(nobsn, size(C,1));
%     x(1,:) = [0 ; par.mu / par.gamma]'; % Initial value of x.
    x(1,:) = [0; par.mu]';

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

vsig = par.s;
% Covariance matrix of measurement error
if correlation == 0
    V = diag(vsig.^2);
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
    D = diag(vsig.^2);
    V = D^(1/2) * CorMat * D^(1/2);

else
    error('correlation must be 0 or 1.')
end

% Distribution of the state error.
normmu = repelem(0, size(x, 2));
normsigma = [1, par.rho_chixi; par.rho_chixi, 1];
burn_in = 500;

% Initial values of state variables.
chi(1) = x(1,1);
xi(1) = x(1,2);

% Simulating state variables
for i = 2:nobsn + burn_in + 1
    z(i-1, :) = mvnrnd(normmu, normsigma);
    chi(i) = chi(i-1) * exp(-par.kappa * deltat) + par.sigmachi * sqrt((1 - exp(-2 * par.kappa * deltat)) / (2 * par.kappa)) * z(i-1, 1);
    if LT == "GBM"
        xi(i) = xi(i-1) + par.mu * deltat + par.sigmaxi * sqrt(deltat) * z(i-1, 2);
    elseif LT == "OU"
        xi(i) = xi(i-1) * exp(-par.gamma * deltat) + par.mu / par.gamma * (1 - exp(-par.gamma * deltat)) + ...
            par.sigmaxi * sqrt((1 - exp(-2 * par.gamma * deltat)) / (2 * par.gamma)) * z(i-1, 2);
    end
end

x = [chi;xi]';
x = x(burn_in + 2:end, :); % Removing burn-in period.

% Simulating measurable variables
vnormmu = repelem(0, ncontracts)';
vnormsigma = V;
y = zeros(nobsn, ncontracts);

% Simulate AR Errors
if n_lags > 0
    AR = {};
    for i = 1:n_lags
        AR{i} = [diag(par.phi(((i-1)*ncontracts+1):ncontracts*i))];
    end
    V_temp = triu(V);
    V_temp = V_temp + V_temp' - diag(diag(V_temp));
    error_mdl = varm('Constant', vnormmu, 'AR', AR, 'Covariance', V_temp);
    v = simulate(error_mdl, nobsn);
elseif n_lags == 0
    v = mvnrnd(vnormmu, vnormsigma, nobsn);
end


for i = 1:nobsn
    y(i,:) = d(:,i) + B(:,:,i) * x(i,:)' + v(i,:)';
end

simY = y;

% ------------------------ End of simulation ----------------------------%