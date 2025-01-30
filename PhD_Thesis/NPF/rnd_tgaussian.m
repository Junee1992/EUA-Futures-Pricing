function X = rnd_tgaussian(m,s2,t0,t1)
%
% function X = rnd_tgaussian(m,s2,t0,t1)
%
% Sampling from a Gaussian with mean m and variance s2, truncated outside the interval [t0,t1], by rejection sampling.
%
% 'm' and 's2' are vectors of the same size, containing the means and variances, respectively.
% 't0' and 't1' are the bounds; they can be either vectors of the same size as 'm' and 's2' or scalars.
%

% First step
X = m + sqrt(s2).*randn(size(m));

% Check for rejection
Idx = find( (X<t0) | (X>t1) );

% Keeps drawing until Idx is empty (no sample is rejected)
while not(isempty(Idx))
	X(Idx) = m(Idx) + sqrt(s2(Idx)).*randn(size(Idx));
	Idx = find( (X<t0) | (X>t1) );
end; %while