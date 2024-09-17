function  alpha = TheilSen(X, y,deltax,deltay);
% THEILSEN performs Theil-Sen robust simple linear regression(s) of X on y.
% 
% For convenience, the input X may contain more than one predictor variable;
% but note that two or more predictor variables in X are treated as independent
% simple regressions: Do not confuse the output with multiple regression.
%
% With x denoting a single column vector of X, the returned slope estimator is
% the median of the set of slopes of data point pairs,
% (y(j) - y(i)) / (x(j) - x(i)), with x(i) ≠ x(j).
% The returned intercept is the median of y - slope_estimator * x.
%
% THEILSEN treats NaNs in X or y as missing values and ignores them.
%
% INPUT
%   X: One or more column vectors containing explanatory/predictor variables.
%      Rows represent observations (assumed to be i.i.d.).
%   y: A column vector containing the observations of the response variable.
%
% OUTPUT
%    coef: Estimated regression coefficients for each predictor column in X,
%          with respect to the response variable y. Each column in coef
%          corresponds to one predictor in X, i.e. it has as many columns as X.
%          The first row, i.e. coef(1, :), contains the estimated offset(s).
%          The second row, i.e. coef(2, :), contains the estimated slope(s).
%          (This output format was chosen to avoid confusion, e.g. with previous
%           versions of this code.)
%   rsqrd: Ordinary R² (coefficient of determination) per predictor column in X.
%
% EXAMPLE
%   See accompanying file example.m.
%
% REFERENCE
%   - Gilbert, Richard O. (1987), "6.5 Sen's Nonparametric Estimator of Slope",
%     Statistical Methods for Environmental Pollution Monitoring,
%     John Wiley and Sons, pp. 217-219, ISBN 978-0-471-28878-7
%   - Sen, P. K. (1968). Estimates of the Regression Coefficient Based on
%     Kendall’s Tau. Journal of the American Statistical Association, 63(324),
%     1379–1389. https://doi.org/10.1080/01621459.1968.10480934
%
% AUTHORS
%   2014-15 Zachary Danziger
%   2022 Johannes Keyser
%
% LICENSE
%   BSD 2-clause "simplified" license, see accompanying file license.txt.

sizeX = size(X);
sizeY = size(y);

if length(sizeY) ~= 2 || sizeY(1) < 2 || sizeY(2) ~= 1 || ~isnumeric(y)
    error('Input y must be a column array of at least 2 observed responses.')
end

if length(sizeX) ~= 2 || ~isnumeric(X)
    error('Input X must be one or more column arrays of predictor variables.')
end

if sizeX(1) ~= sizeY(1)
    error('The number of rows (observations) of X and y must match.')
end

Num_Obs = sizeX(1);  % rows in X (and y) are observations
Num_Pred = sizeX(2);  % columns in X are (independent) predictor variables


% calculate slope for all pairs of data points
C = nan(Num_Obs, Num_Obs);
for i = 1:Num_Obs-1
      dy = y(i) - y(i:end);
      dx = (X(i) - X(i:end));
      C(i, i:end) = dy ./ dx;
      C(i, find(((abs(dx)>deltax).*(abs(dy)>deltay))==0)) = NaN;

end
% relabel infinite values (due to div by 0 = x(i) - x(j)) as NaNs-to-ignore.
C(isinf(C)) = NaN;
% estimate slope as the median of all pairwise slopes
b1 = median(C(:), 'omitnan');
% estimate offset as the median of all pairwise offsets
b0 = median(y - b1 * X, 'omitnan');
alpha = [ b1 ; b0 ]

end
