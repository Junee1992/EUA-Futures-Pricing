%% Sets parameters as a structure

function param = set_parameters(LT, ncontracts, par_names, par, correlation, serial)
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

% if startsWith(seasonal, 'sinusoid_')
%     j = 1;
%     for i = [find(startsWith(par_names, 'a_'))]
%         param.a(j) = par_init(i);
%         j = j+1;
%     end
% 
%     j = 1;
%     for i = [find(startsWith(par_names, 'b_'))]
%         param.b(j) = par_init(i);
%         j = j+1;
%     end
% 
%     if seasonal == "sinusoid_b"
%         j = 1;
%         for i = [find(startsWith(par_names, 'A_'))]
%             param.A(j) = par_init(i);
%             j = j+1;
%         end
%     end
% end

if serial == "yes"
    j = 1;
    for i = min(find(startsWith(par_names, 'phi_'))):max(find(startsWith(par_names, 'phi_')))
        param.phi(j) = par(i);
        j = j+1;
    end
elseif serial == "no"
    param = param;
end



