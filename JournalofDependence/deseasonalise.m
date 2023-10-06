function [season, freq, fitresult, gof] = deseasonalise(price, n_season, seasonality)

y = mean(price,2);
nobsn = size(y,1);
x = [1:nobsn]';

[pxx, w] = pwelch(y);
pxx1 = 10*log10(pxx);
w1 = w/pi;

% pxx1_temp = pxx1;
% periodicity = [];
% pxx_max = [];
% 
% if n_season == 0
%     periodicity = [];
% elseif n_season > 0
%     pxx_max(1) = max(pxx1_temp);
%     index(1) = find(pxx1_temp == max(pxx1_temp));
%     periodicity(1) = w1(index(1));
%     pxx1_temp(1:index(1)) = NaN;
% end
% 
% if n_season > 1
%     for i = 2:n_season
%         pxx_max(i) = max(pxx1_temp(pxx1_temp < pxx_max(i-1)));
%         index(i) = find(pxx1_temp == pxx_max(i));
%         periodicity(i) = w1(index(i));
%         pxx1_temp(1:index(i)) = NaN;
%     end
% end

[pks, freq] = findpeaks(pxx1, w1);

a = [];
b = [];

for i = 1:n_season
    a = [a sprintf("a%d", i)];
    b = [b sprintf("b%d", i)];
end

txt = [];

if seasonality == "sinusoid_a"
    txt = [sprintf(a{1}), sprintf('*cos(1*2*pi*x)+'), sprintf(b{1}), sprintf('*sin(1*2*pi*x)')];
    for i = 2:n_season
        txt = [txt, sprintf('+'), sprintf(a{i}), sprintf('*cos('), sprintf(num2str(i)), sprintf('*2*pi*x)+'), sprintf(b{i}), sprintf('*sin('), sprintf(num2str(i)), sprintf('*2*pi*x)')];
    end
    ft = fittype(txt, 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';


elseif seasonality == "sinusoid_b"
    txt = [sprintf(a{1}), sprintf('*cos('), sprintf(num2str(freq(1))), sprintf('*2*pi*x)+'), sprintf(b{1}), sprintf('*sin('), sprintf(num2str(freq(1))), sprintf('*2*pi*x)')];
    for i = 2:n_season
        txt = [txt, sprintf('+'),sprintf(a{i}), sprintf('*cos('), sprintf(num2str(freq(i))), sprintf('*2*pi*x)+'), sprintf(b{i}), sprintf('*sin('), sprintf(num2str(freq(i))), sprintf('*2*pi*x)')];
    end
    ft = fittype(txt, 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.StartPoint = [repelem(0, n_season*2)];
    opts.Display = 'Off';

end

% if seasonality == "sine"
% %     txt = 'a1*sin(b1*x+c1)';
% %     opts.Lower = [-Inf 0 -Inf];
% %     opts.StartPoint = [0.349994354999 0.00587762891223535 -1.40544117299823];
% %     ft = fittype( 'sin1' );
% %     opts = fitoptions('Method', 'NonlinearLeastSquares');
% %     opts.Display = 'Off';
% 
%     ft = fittype( 'fourier1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [0 0 0 0.00587762891223535];

% end


[fitresult, gof] = fit(x, y, ft, opts );

season = fitresult(x);
y_deseason = y - season;

% figure;
% plot(y);
% hold on
% plot(seasonality);
% plot(y_deseason);
% legend('Mean of Futures Price', 'Seasonality', 'Deseasonalised')
