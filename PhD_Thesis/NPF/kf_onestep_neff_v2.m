function [lW, particle_post, particle, wu, w, Neff]  = kf_onestep_neff_v2(y, ttm, par_particle, particle, model_options, part_est, num)
% [lW, particle_post, particle] = kf_onestep(y, ttm, par_particle,
% particle, model_options)
%
% Performs SIR at a time point t.
%
% Inputs:
% y = an [n x K] matrix, n = number of observations, K = number of
% contracts
% ttm = an [n x K] matrix representing the time to maturity
% par_particle = a vector of parameters in a single external particle
% particle = an [N x 2] matrix representing inner particles
% model_options = options for the models
%
% Outputs: For m-th particle, m = 1, ..., M,
% lW = mean of "N" log-likelihoods to represent the likelihood for
% particle_post = weighted average of particles
% particle = resampled particles to be used in the next time t+1.

% rng(1000);
fields = fieldnames(model_options);
for i = 1:numel(fields)
    fieldName = fields{i};
    eval([fieldName ' = model_options.(fieldName);']);
end

ncontracts = size(y,2);
par = set_parameters(LT, ncontracts, par_names, par_particle, correlation, "no", err);

C = [0 ; (par.mu / par.gamma) * (1 - exp(-par.gamma * deltat))];
G = [exp(-par.kappa * deltat), 0; 0, exp(-par.gamma * deltat)];
w11 = (1 - exp(-2 * par.kappa * deltat)) / (2 * par.kappa) * par.sigmachi^2;
w12 = (1 - exp(-(par.kappa + par.gamma) * deltat)) / (par.kappa + par.gamma) * (par.sigmachi * par.sigmaxi * par.rho_chixi);
w22 = (1 - exp(-2 * par.gamma * deltat)) / (2 * par.gamma) * par.sigmaxi^2;
W = [w11, w12; w12, w22];

[ev, eva] = eig(W);
eva(eva < 0) = 0;
W = ev * eva * ev';


d1 = (1 - exp(-2 * par.kappa * ttm)) * par.sigmachi^2 / (2 * par.kappa);
d2 = (1 - exp(-2 * par.gamma * ttm)) * par.sigmaxi^2 / (2 * par.gamma);
d3 = (1 - exp(-(par.kappa + par.gamma) * ttm)) * 2 * par.sigmachi * par.sigmaxi * par.rho_chixi / (par.kappa + par.gamma);
d = (par.mu - par.lambdaxi) / par.gamma * (1 - exp(-par.gamma * ttm)) - (par.lambdachi / par.kappa) * (1 - exp(-par.kappa * ttm)) + (1/2) * (d1 + d2 + d3);
d = d';

B1 = exp(-par.kappa * ttm);
B2 = exp(-par.gamma * ttm);
B = [B1; B2]';

vsig2 = par.s.^2;
if correlation == 0
    V = diag(vsig2);
elseif correlation == 1
    correl = par.rho;

    % Manually creating the correlation matrix of measurement errors
    CorMat = diag(repelem(1, ncontracts));
    for a = 1:ncontracts
        for b = 1:ncontracts
            if a == b
                CorMat(a,b) = 1;
            else
                CorMat(a,b) = correl(a) * correl(b);
            end
        end
    end
    D = diag(vsig2);
    V = D^(1/2) * CorMat * D^(1/2);
end

    
% Likelihood
if err == "normal"
    parfor j = 1:N
        % State Prediction
        particle(j,:) = (C + G * particle(j,:)' + mvnrnd(repelem(0, size(C,1)),W)')';
        % particle(j,:) = (C + G*particle(j,:)')';
        yfit(j,:) = d + B * particle(j,:)';
        z_esti(:,j) = d + B * particle(j,:)' + mvnrnd(repelem(0, ncontracts), V)';
        ett(j,:) = y - z_esti(:,j)';
        % lkhd(j) = log(likelihood(y, z_esti(:,j)', V));
        lkhd(j) = log(likelihood(ett(j,:), repelem(0, ncontracts), V));
    end
elseif err == "laplace"
    parfor j = 1:N
        % State Prediction
        particle(j,:) = (C + G * particle(j,:)' + mvnrnd(repelem(0, size(C,1)),W)')';
        yfit(j,:) = d + B * particle(j,:)';
        z_esti(j,:) = d + B * particle(j,:)' + randLaplace(repelem(0, ncontracts), V, par.sG, 1)';
        ett(j,:) = y - z_esti(j,:);
        lkhd(j) = log(mvLaplacedensity(ett(j,:)', [repelem(0, ncontracts)]', B*W*B'+V, par.sG));
    end
elseif err == "hyperbolic"
    parfor j = 1:N
        % State Prediction
        particle(j,:) = (C + G * particle(j,:)' + mvnrnd(repelem(0, size(C,1)),W)')';

        z_esti = d + B * particle(j,:)' + randmgh(par.lambda, par.psi, par.chi, repelem(0, ncontracts)', V, par.nu');
        ett = y - z_esti';
        lkhd(j) = log(dghypmv(ett', par.lambda, par.psi, par.chi, repelem(0, ncontracts)', V, par.nu'));
    end
end


% Sample impoverishment
if not( isempty( find( isnan(lkhd) | isinf(lkhd), 1 ) ) )
    idx_bad = find( isnan(lkhd) | isinf(lkhd) );
    lkhd( idx_bad ) = -1000;
    % lkhd( idx_bad ) = min(lkhd(~isinf(lkhd) & ~isnan(lkhd))) - 100;
    v1 = part_est(:,1); v2 = part_est(:,2);
    % particle(idx_bad,:) = [mean(part_est(:,1)) * ones(size(idx_bad,1), 1), mean(part_est(:,2)) * ones(size(idx_bad,1), 1)] +...
    %     0.1 * randn(size(idx_bad',1), 2);
    particle(idx_bad,:) = [mean(v1(v1 ~= 0)) * ones(size(idx_bad,1), 1), mean(v2(v2 ~= 0)) * ones(size(idx_bad,1), 1)] +...
         0.01 * randn(size(idx_bad',1), 2);

% else
%     pause(10);
end

% Calculate the weights
wu = exp(lkhd - max(lkhd));
ind = wu < realmin;
wu(ind) = realmin;
% wu = lkhd;
w = wu/sum(wu); % Normalise
Neff = 1/sum(w.^2);
% 
% for j = 1:7
% subplot(4,2,j)
% yline(y(:,j));
% hold on;
% % plot(repelem(num, N), yfit(:,j),'r.');
% % errorbar(num, mean(z_esti(j,:)'), mean(z_esti(j,:)')-min(z_esti(j,:)'), max(z_esti(j,:)')-mean(z_esti(j,:)'));
% % errorbar(num, mean(yfit(j,:)'), mean(yfit(j,:)')-min(yfit(j,:)'), max(yfit(j,:)')-mean(yfit(j,:)'));
% [f, xi] = ksdensity(z_esti(j,:)', 'Function', 'pdf');
% plot(xi, f);
% hold off;
% end
% % subplot(4,2,8)
% % hold on
% % plot(num, Neff, '*')
% % hold on
% pause(0.01);
% 
% % % if Neff < 0.5*N
% % %     % Resampling
% % %     [particle_resamp, ~] = Resample(particle', wu);
% % %     particle = particle_resamp';
% % % end

if Neff<0.3*N
    % try
     [mu, si, p] = EMGM(particle', w', 5);
    % catch ME
    %  save(sprintf('debug_iter.mat'), 'particle', 'w', 'par', 'C', 'G', 'W', 'd', 'B');
    % end

     nw=length(p);
     % nw_i=[nw_i,nw]; 
     for j=1:N
         indw  = randsample(nw, 1, true,p);
         particle(j,:) = mvnrnd(mu(:,indw), si(:,:,indw));
     end
     % w=1/N*ones(N,1)'; % reweighting
%     % H2=-w'*log(mvnpdf(part_u,mu_u,cov_u)+realmin);
%     % H=[H;H1,H2];
% 
 end
% if Neff>0.3*N
%     "hi";
% end
lW = log(mean(wu)) + max(lkhd); % Average of weights representing the likelihood for the m-th external particle
% lW = log(mean(wu));
% lW = mean(wu);
particle_post = w * particle; % Weighted average of x_t.

% Resampling
%[particle_resamp, ~] = Resample(particle', wu);
%particle = particle_resamp';
end