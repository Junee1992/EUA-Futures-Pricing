function [varargout] = crps(fcst,obs,varargin)
% CRPS  Calculates continuous rank probability score (CRPS)
% 
% The CRPS measures the closeness of forecast distribution (fcst) and
% corresponding observation (obs). For more info see
% Jolliffe, I. T. and Stephenson, D. B. (eds.) 2012. Forecast verification:
% a practitioner's guide in atmospheric science, West Sussex, England: Wiley.
%
% USAGE:
%        [mean_CRPS] = crps(fcst,obs);
%        [mean_CRPS] = crps(fcst,obs,plot_pos);
%        [mean_CRPS,crps_values,num] = crps(fcst,obs);
%
% INPUT:
%   obs        -   Vector of observations
%   fcst       -   Matrix of Ensemble forecast of size N x M. NB: N
%                  must equal length(obs), M equals the number of
%                  ensemble members.
%  plot_pos    -   plotting positions that determine cumulative distribution function.
%                  Options are as follows:
%  Param        Equation            Description
%  'ecdf'       i/n                 linear interpolation (emprical cdf e.g. ecdf in
%                                   MATLAB, Default)
%  'Weibull'    i/(n+1)             Unbiased exceedance probabilities for
%                                   all distributions
%  'Median'     (i-0.3175)/(n+0.365)Median exceedance probabilities for
%                                   all distributions
%  'APL'        (i-0.35)/n          Used with PWMs.
%  'Blom'       (i-0.375)/(n+0.25)  Unbiased normal quantiles
%  'Cunnane'    (i-0.4)/(n+0.2)     Approx. quantile-unbiased
%  'Gringorten' (i-0.44)/(n+0.12)   Optimised for Gumbel distribution
%  'Hazen'      (i-0.5)/n           [Default] a traditional choice
%  'Bernad'     (i-0.3)/(n + 0.4)   Benard, A and Bos-Levenbach, E. C. authors, 1953
%  [p1 p2]      (i-p1)/(n+p2)       2-element vector defining custom values
%                                   for p1 & p2. p1 & p2 must take values
%                                   between 0 and 1.
% See: Stedinger JR, Vogel RM, Foufoula-Georgiu E (1995) 'Frequency
% Analysis of Extreme Events' in Maidment, D (ed.) Handbook of Hydrology.
% New York: McGraw-Hill.
%
% OUTPUT:%   
%   mean_CRPS     -  Mean of non missing CRPS values
%   crps_values    - A vector (length n) of CRPS values
%   num      -       number of non missing CRPS values used to compute mean_CRPS
%
% EXAMPLES:
%    fcst = rand(1000,1000);
%    obs = rand(1000,1);
%    [meanCRPS] = crps(fcst,obs);
%    [meanCRPS crps_values] = crps(fcst,obs,'Bernad');
%
% See also:
% Author: Durga Lal Shrestha
% CSIRO Land & Water, Highett, Australia
% eMail: durgalal.shrestha@gmail.com
% Website: https://sites.google.com/site/durgalalshrestha/
% Copyright 2014 Durga Lal Shrestha
% $First created: 10-Sep-2014
% $Revision: 1.0.0 $ $Date: 10-Sep-2014 14:50:33 $
% ***********************************************************************
%% INPUT ARGUMENTS CHECK
narginchk(2,3)
plot_pos = 'ecdf';
if nargin>2
    plot_pos = varargin{1};
end
% Check ploting position method
if ischar(plot_pos)
    if strcmpi(plot_pos,'ecdf')
        p1 = 0; p2 = 0;
    elseif strcmpi(plot_pos,'Weibull')
        p1 = 0; p2 = 1;
    elseif strcmpi(plot_pos,'Median')
        p1 = 0.3175; p2 = 0.365;
    elseif strcmpi(plot_pos,'APL')
        p1 = 0.35; p2 = 0;
    elseif strcmpi(plot_pos,'Blom')
        p1 = 0.375; p2 = 0.25;
    elseif strcmpi(plot_pos,'Cunnane')
        p1 = 0.4; p2 = 0.2;
    elseif strcmpi(plot_pos,'Gringorten')
        p1 = 0.44; p2 = 0.12;
    elseif strcmpi(plot_pos,'Hazen')
        p1 = 0.5; p2 = 0;
    elseif strcmpi(plot_pos,'Bernad')
        p1 = 0.3; p2 = 0.4;
    else
        error('CRPS:BadPlottingPosition','Plotting position not recognised');
    end
else
    if length(plot_pos)~=2
        error('CRPS:BadPlottingPosition','Plotting position vector does not have 2 elements');
    end
    if max(plot_pos)>1 || min(plot_pos)<0
        error('CRPS:BadPlottingPosition','Plotting position parameters must be between one and zero');
    end
    p1 = plot_pos(1); p2 = plot_pos(2);
end
if ~ismatrix(fcst) || ~isvector(obs) 
    error('CRPS:BadDimension','CRPS works only for two dimensional forecast matrix and vector (or scalar) observation');
end
% Figure out which dimension crps will work along.
sz_fcast = size(fcst);
m = find(length(obs)==sz_fcast);
if m==2
    fcst=fcst';
end
obs = obs(:);
[n, m]=size(fcst);
if n~=length(obs)
    error('CRPS:BadDimension','number of forecasts should be equal to the length of observation')
end
% Treat negatives are missing values
fcst(fcst<0)=NaN;
obs(obs<0)=NaN;
% Compute probabilty
pi = ((1:m)-p1)./(m+p2);
p2i = pi.^2;  %ith probability squared
minusP2i = (1.0-pi).^2; %
crps_values = NaN(n,1);
fcst = sort(fcst,2);
% Loop through all number of observations
for i=1:n
    % check if there are any missing values
    missingFcst = any(isnan(fcst(i,:)));
    missingObs = isnan(obs(i));
    if  ~missingFcst && ~missingObs
        ind = find(fcst(i,:)<=obs(i),1,'last');
        if ~isempty(ind) % i.e. obs(i) >  fcst(i,1)
            % left of the observation
            crpsLeft = 0;
            if ind>1
                fcstLeft = fcst(i,1:ind);
                dxLeft = diff(fcstLeft);
                pLeft = p2i(1:ind-1);
                crpsLeft = pLeft*dxLeft';
            end
            
            if obs(i)<fcst(i,end)
                % right of the observation
                fcstRight = fcst(i,ind+1:end);
                dxRight = diff(fcstRight);
                if isempty(dxRight)
                    crpsRight = 0;
                else
                    pRight = minusP2i(ind+1:end-1);
                    crpsRight = pRight*dxRight';
                end
                % when the cdf crosses the observation
                % left part
                crpsCentreLeft = p2i(ind).*(obs(i)-fcst(i,ind));
                % right part
                crpsCentreRight = minusP2i(ind).*(fcst(i,ind+1)-obs(i));
                crps_values(i) = crpsLeft + crpsRight+crpsCentreLeft+crpsCentreRight;
                
            else  % if observation lies right of all members of forecast (ie. obs>fcst(end))
                crps_right_outside = 1.0^2*(obs(i)-fcst(i,end));
                crps_values(i)  =  crps_right_outside + crpsLeft;
            end
            
        else  % observation lies left of the all member of forecast (ie. obs(i) < fcst(i,1))
            dxRight = diff(fcst(i,:));
            pRight = minusP2i(1:end-1);
            crpsRight = pRight*dxRight';
            crps_left_outside = 1.0^2*(fcst(i,1)-obs(i));
            crps_values(i)  =  crps_left_outside+crpsRight;
        end
    end  % not missing
end  % all forecasts
%% Output
if n==1
    varargout{1} = crps_values;
else
    ind = isnan(crps_values);
    cpNotNaN = crps_values(~ind);
    varargout{1} = mean(cpNotNaN);
end
if nargout>1    
    varargout{2}= crps_values;
    varargout{3}= length(cpNotNaN);
end
