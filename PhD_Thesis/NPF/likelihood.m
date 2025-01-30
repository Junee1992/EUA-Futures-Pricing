%% The Multivariate-Normal likelihood function
% function p = likelihood(Z_est, Z_mea, V)
% Inputs:
% Z_est : an [n x k] matrix representing the estimates of Z
% Z_mea : an [n x k] matrix representing the actual values of Z
% V : a [k x k] covariance matrix.
%
% Output:
% p : the likelihood

function [ p ] = likelihood( Z_est, Z_mea, V)
    k = size(Z_est,2);
    Z_mea = Z_mea'; Z_est = Z_est';
    p = (2*pi)^(-k/2) * det(V)^(-1/2) * exp(-1/2 * (Z_mea - Z_est)' * inv(V) * (Z_mea - Z_est));
end

