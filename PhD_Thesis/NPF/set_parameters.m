%% Sets parameters as a structure

function param = set_parameters(LT, ncontracts, par_names, par, correlation, serial, err)
par_names = par_names';
param = [];

if LT  == "GBM"
    for i = 1:7
        param.(genvarname(par_names(i))) = par(i);
    end

    j = 1;
    for i = 8:7+ncontracts
        param.s(j) = par(i);
        j = j+1;
    end

    if correlation == 0
        param = param;

    elseif correlation == 1
        j = 1;
        for i = 8+ncontracts : 7+2*ncontracts
            param.rho(j) = par(i);
            j = j+1;
        end

    else
        error('correlation must be 0 or 1.')
    end

elseif LT == "OU"
    for i = 1:8
        param.(genvarname(par_names(i))) = par(i);
    end

    j = 1;
    for i = 9:8+ncontracts
        param.s(j) = par(i);
        j = j+1;
    end

    if correlation == 0
        param = param;
        param.rho = [];

    elseif correlation == 1
        j = 1;
        for i = 9+ncontracts : 8+2*ncontracts
            param.rho(j) = par(i);
            j = j+1;
        end

    else
        error('correlation must be 0 or 1.')
    end
end

if serial == "yes"
    j = 1;
    for i = min(find(startsWith(par_names, 'phi_'))):max(find(startsWith(par_names, 'phi_')))
        param.phi(j) = par(i);
        j = j+1;
    end
elseif serial == "no"
    param = param;
    param.phi = [];
end

if err == "normal"
    param = param;
elseif err == "laplace"
    param.sG = par(end);
elseif err == "hyperbolic"
    param.lambda = par(8+length(param.s)+length(param.rho) + length(param.phi) + 1);
    param.psi= par(8+length(param.s)+length(param.rho) + length(param.phi) + 2);
    param.chi = par(8+length(param.s)+length(param.rho) + length(param.phi) + 3);
    param.nu = par(8+length(param.s) + length(param.rho) + length(param.phi) + 4:8+length(param.s) + length(param.rho) + length(param.phi) + 3+ncontracts);
end


