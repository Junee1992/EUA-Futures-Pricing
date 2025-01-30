function negLogLikelihood = GHobjectiveFunction(theta, ett, cov_matrix)
    [nobsn, ncontracts] = size(ett);
    % Extract alpha and the diagonal elements (sigma9, ..., sigma13)
    lambda = theta(1); psi = theta(2); chi = theta(3);
    nu = theta(4:end);
    % if length(theta(1:end-1)) > ncontracts
    %     correl = theta(ncontracts+1:2*ncontracts);
    % 
    %     % Manually creating the correlation matrix of measurement errors
    %     CorMat = diag(repelem(1, ncontracts));
    %     for a = 1:ncontracts
    %         for b = 1:ncontracts
    %             if a == b
    %                 CorMat(a,b) = 1;
    %             else
    %                 CorMat(a,b) = correl(a) * correl(b);
    %             end
    %         end
    %     end
    %     D = diag(theta(1:ncontracts).^2);
    %     cov_matrix = D^(1/2) * CorMat * D^(1/2);
    % else
    %     correlation = 0;
    %     cov_matrix = diag(theta(1:end-1).^2);
    % end

    lkhd = 0;
    for i = 1:nobsn
        % Create the diagonal covariance matrix
        
       
        % Compute the log-likelihood based on the Multivariate Laplace density
        lkhd = lkhd + log(dghypmv(ett(i,:)', lambda, psi, chi, repelem(0, ncontracts)', cov_matrix, nu'));
    end
    negLogLikelihood = -lkhd;  % Return the negative of the log-likelihood
end
