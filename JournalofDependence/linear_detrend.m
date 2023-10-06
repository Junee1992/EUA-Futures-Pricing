function [y_detrend, trend, fitresult] = linear_detrend(y)

[nobsn, ncontracts] = size(y);
t = [1:nobsn]';
tt = repmat(t, ncontracts, 1);
yy = reshape(y, ncontracts*nobsn, 1);

ft = fittype('alpha+beta*x', 'independent', 'x', 'dependent', 'y');
opts = fitoptions('Method', 'NonLinearLeastSquares');
opts.StartPoint = [1, 1];
[fitresult, gof] = fit(tt, yy, ft, opts );

y_detrend = y - fitresult(t);
trend = fitresult(t);