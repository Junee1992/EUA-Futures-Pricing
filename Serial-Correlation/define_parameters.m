% Define Parameters

function par_names = define_parameters(LT, ncontracts, correlation, n_lags, test)

par_names= [];

if LT == "GBM"
    if test == "yes"
        par_names = ["kappa", "sigmachi", "lambdachi", "mu", "sigmaxi", "rnmu", "rho_chixi"];
    elseif test == "no"
        par_names = ["kappa", "sigmachi", "lambdachi", "mu", "sigmaxi", "lambdaxi", "rho_chixi"];
    end
elseif LT == "OU"
    par_names = ["kappa", "sigmachi", "lambdachi", "gamma", "mu", "sigmaxi", "lambdaxi", "rho_chixi"];
else
    error('Please specify the process of the long-term factor LT.')
end

for i = 1:ncontracts
    par_names = [par_names, sprintf("s_%d", i)];
end

if correlation == 0
    par_names = par_names;

elseif correlation == 1
    for i = 1:ncontracts
        par_names = [par_names, sprintf("rho_%d", i)];
    end

else
    error('correlation must be 0 or 1.')
end

if n_lags == 0
    par_names = par_names;
elseif n_lags > 0
    for j = 1:n_lags
        for i = 1:ncontracts
            par_names = [par_names, sprintf("phi_%d%d", j, i)];
        end
    end
end

% if seasonal == "sinusoid_a"
%     if n_season > 0
%         par_names = [par_names, "a_1", "b_1"];
%     end
%     if n_season > 1
%         for i = 2:n_season
%             par_names = [par_names, sprintf("a_%d", i), sprintf("b_%d", i)];
%         end
%     end
% 
% elseif seasonal == "sinusoid_b"
%     if n_season > 0
%         par_names = [par_names, "a_1", "b_1", "A_1"];
%     end
%     if n_season > 1
%         for i = 2:n_season
%             par_names = [par_names, sprintf("a_%d", i), sprintf("b_%d", i), sprintf("A_%d", i)];
%         end
%     end
% 
% elseif seasonal == "no"
%     par_names = par_names;
% 
% else
%     error('Please select available function for seasonality.')
% end



