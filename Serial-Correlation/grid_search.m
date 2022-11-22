function init = grid_search(lb, ub, ngrid, npar_grid, n_par, n_lag, ncontracts, LT)

% Grid search
mid = (lb + ub) / 2;
mp1 = (mid + lb) / 2;
mp2 = (mid + ub) / 2;
grid = [mp1; mid; mp2]';

% est = zeros(ngrid^npar_grid, n_par+1);
init = [];

if LT == "GBM"
    for k = grid(1,:)
        for sc = grid(2,:)
            for m = grid(4,:)
                for sx = grid(5,:)
                    for rh = grid(6,:)
                        init = [init; k, sc, 0, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.8, n_lag)];
                    end
                end
            end
        end
    end

elseif LT == "OU"
    for k = grid(1,:)
        for sc = grid(2,:)
            for g = grid(4,:)
                for m = grid(5,:)
                    for sx = grid(6,:)
                        for rh = grid(8,:)
                            init = [init; k, sc, 0, g, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.8, n_lag)];
                        end
                    end
                end
            end
        end
    end
end