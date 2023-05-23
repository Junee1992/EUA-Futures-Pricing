%% Constraints of model parameters

function [A, b, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, n_lag)

A = zeros(1, n_par);
A(1) = -1; A(4) = 1;
b = 0;
lb = [];
ub = [];

if LT == "GBM"
    lb(1:7) = [1e-5, 1e-5, -10, -10, 1e-5, -10, -1];
    lb(8:7+ncontracts) = 1e-3;
    ub(1:7) = [10, 10, 10, 10, 10, 10, 1];
    ub(8:7+ncontracts) = 1;
    if correlation == 1
        lb(8+ncontracts:7+2*ncontracts) = 0;
        ub(8+ncontracts:7+2*ncontracts) = 1;
    end

elseif LT == "OU"
    lb(1:8) = [1e-5, 1e-5, -5, 1e-5, -5, 1e-5, -5, -1];
    lb(9:8+ncontracts) = 1e-3;
    ub(1:8) = [5, 1, 1, 5, 10, 1, 1, 1];
    ub(9:8+ncontracts) = 1;
    if correlation == 1
        lb(9+ncontracts:8+2*ncontracts) = 0;
        ub(9+ncontracts:8+2*ncontracts) = 1;
    end
   
end

if n_lag > 0
    lb = [lb, repelem(-0.9999, n_lag*ncontracts)];
    ub = [ub, repelem(0.9999, n_lag*ncontracts)];
end