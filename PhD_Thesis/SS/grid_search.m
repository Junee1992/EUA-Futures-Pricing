function init = grid_search(lb, ub, n_lag, ncontracts, LT, correlation)

% Grid search
mid = (lb + ub) / 2;
mp1 = (mid + lb) / 2;
mp2 = (mid + ub) / 2;
grid = [lb; mp1; mid; mp2; ub]';

% est = zeros(ngrid^npar_grid, n_par+1);
init = [];

if LT == "GBM"
    for k = grid(1,:)
        for sc = grid(2,:)
            for m = grid(4,:)
                for sx = grid(5,:)
                    for rh = grid(7,:)
                        if correlation == 0
                                if n_lag > 1
                                    init = [init; k, sc, 0, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.7, ncontracts), repelem(0.1, ncontracts*(n_lag-1))];
                                elseif n_lag == 1
                                    init = [init; k, sc, 0, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.7, ncontracts)];
                                elseif n_lag == 0
                                    init = [init; k, sc, 0, m, sx, 0, rh, repelem(0.01, ncontracts)];
                                end
                            elseif correlation == 1
                                if n_lag > 1
                                    init = [init; k, sc, 0, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.9, ncontracts), repelem(0.7, ncontracts), repelem(0.1, ncontracts*(n_lag-1))];
                                elseif n_lag == 1
                                    init = [init; k, sc, 0, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.9, ncontracts), repelem(0.7, ncontracts)];
                                elseif n_lag == 0
                                    init = [init; k, sc, 0, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.9, ncontracts)];
                                end
                            end
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
                            if correlation == 0
                                if n_lag > 1
                                    init = [init; k, sc, 0, g, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.7, ncontracts), repelem(0.1, ncontracts*(n_lag-1))];
                                elseif n_lag == 1
                                    init = [init; k, sc, 0, g, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.7, ncontracts)];
                                elseif n_lag == 0
                                    init = [init; k, sc, 0, g, m, sx, 0, rh, repelem(0.01, ncontracts)];
                                end
                            elseif correlation == 1
                                if n_lag > 1
                                    init = [init; k, sc, 0, g, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.7, ncontracts), repelem(0.7, ncontracts), repelem(0.1, ncontracts*(n_lag-1))];
                                elseif n_lag == 1
                                    init = [init; k, sc, 0, g, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.7, ncontracts), repelem(0.7, ncontracts)];
                                elseif n_lag == 0
                                    init = [init; k, sc, 0, g, m, sx, 0, rh, repelem(0.01, ncontracts), repelem(0.7, ncontracts)];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end