%% Constraints of model parameters

function [A, b, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, n_lag)

A = zeros(1, n_par);
if LT == "OU"
    A(1,1) = -1; A(1,4) = 1;
end
b = 0;
lb = [];
ub = [];

if LT == "GBM"
    lb(1:7) = [1e-5, 0, -10, -10, 0, -10, -1];
    lb(8:7+ncontracts) = 1e-10;
    ub(1:7) = [10, 10, 10, 10, 10, 10, 1];
    ub(8:7+ncontracts) = 1;
    if correlation == 1
        lb(8+ncontracts:7+2*ncontracts) = -0.9999;
        ub(8+ncontracts:7+2*ncontracts) = 0.9999;
    end

elseif LT == "OU"
    lb(1:8) = [1e-5, 1e-5, -1, 1e-5, -3, 1e-5, -1, -1];
    lb(9:8+ncontracts) = 1e-3;
    ub(1:8) = [2, 1, 1, 1, 3, 2, 1, 1];
    ub(9:8+ncontracts) = 1;
    % lb(1:8) = [0 0 -Inf 0 -Inf 0 -Inf -1];
    % ub(1:8) = [Inf Inf Inf Inf Inf Inf Inf 1];
    % lb(9:8+ncontracts) = 0.000001; ub(9:8+ncontracts) = Inf;

    if correlation == 1
        lb(9+ncontracts:8+2*ncontracts) = -1+1e-4;
        ub(9+ncontracts:8+2*ncontracts) = 1-1e-4;
    end
   
end

if n_lag > 0
    lb = [lb, repelem(-1.9999, n_lag*ncontracts)];
    ub = [ub, repelem(1.9999, n_lag*ncontracts)];
    % for j = 1:ncontracts
    %     A(j+1,9+2*ncontracts+(n_lag-1)*j-1) = 1;
    %     if n_lag > 1
    %         for p = 2:n_lag
    %             A(j+1, 9+2*ncontracts+(n_lag-1)*j-1+(p-1)*ncontracts) = 1;
    %         end
    %     end
    % end
    % b = [b, repelem(1, ncontracts)]';
end

% if n_lag > 0
%     % Ensure AR coefficients are within the range for stationarity
%     lb = [lb, repelem(-0.9999, n_lag*ncontracts)];
%     ub = [ub, repelem(0.9999, n_lag*ncontracts)];
% 
%     % Additional constraint: sum of absolute values of AR coefficients < 1
%     for j = 1:ncontracts
%         A(j+1, 9+2*ncontracts+(n_lag-1)*j-1) = 1;
%         if n_lag > 1
%             for p = 2:n_lag
%                 A(j+1, 9+2*ncontracts+(n_lag-1)*j-1+(p-1)*ncontracts) = 1;
%             end
%         end
%     end
%     b = [b, repelem(1, ncontracts)]';
% 
%     % Add constraint for sum of absolute values of AR coefficients
%     for j = 1:ncontracts
%         A = [A; zeros(1, size(A, 2))];
%         for p = 1:n_lag
%             A(end, 9+2*ncontracts+(n_lag-1)*j-1+(p-1)*ncontracts) = 1;
%         end
%         b = [b; 1];
%     end
% end