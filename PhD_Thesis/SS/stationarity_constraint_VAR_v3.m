function [c, ceq] = stationarity_constraint_VAR_v3(params, k, orders)
    % Nonlinear constraint: Ensure stationarity of each marginal AR(p) process
    % params: vector of AR coefficients
    % k: number of variables (processes)
    % orders: vector of AR orders for each process

    % Initialize indices
    start_idx = 1;
    c = [];
    ceq = []; % No equality constraints

    % Loop through each variable
    for i = 1:k
        p = orders(i); % AR order for the i-th process

        % Extract AR coefficients for the i-th variable
        phi = params(start_idx:(start_idx + p - 1));
        start_idx = start_idx + p;

        % Pad phi with zeros if necessary (for AR(1) treated as AR(2))
        if length(phi) < max(orders)
            phi = [phi, zeros(1, max(orders) - length(phi))];
        end

        % Characteristic polynomial for this marginal process
        characteristic_poly = fliplr([1, -phi]); % 1 - phi1*z - phi2*z^2 - ... - phip*z^p

        % Compute roots of the characteristic equation
        roots_poly = roots(characteristic_poly);

        % Add stationarity constraint for this variable
        c = [c; 1 - abs(roots_poly)]; % 1 - |z| <= 0 ensures |z| > 1.
    end
end
