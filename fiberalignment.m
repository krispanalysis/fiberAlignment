%Written by Kris Pan
%Date: Nov 24th, 2018
%Reference: J Biomed Mater Res A. 2016 Jul;104(7):1680-6. doi: 10.1002/jbm.a.35697. Epub 2016 Mar 17.
%Description: detect edges in an image and create an polar plot of
%occurence of the angles


% enter filename as string (e.g., image.tif);
fileName = input('Enter the name of the file: \n','s')
% read the image
im = imread(fileName);
% crop the image to 1024 x 690. This value should be changed
% according the the size of the image and the part that wants
% to be cropped out. For details refer to
% https://www.mathworks.com/help/images/ref/imcrop.html?searchHighlight=imcrop&s_tid=doc_srchtitle
im = imcrop(im,[0 0 1024 690]);   
% use sobel method to emphasize the edge of the image
filter = fspecial('sobel');
% apply filter 
im = imfilter(im,filter,'replicate');
% get two directional gradients, gradientX and gradientY
[gradientX,gradientY] = imgradientxy(im);
% enter y to plot gradient images
gradientsPlot = input('Would you like to plot the two directional gradients? [y]\n','s');
if gradientsPlot == 'y';
    imshowpair(gradientX, gradientY, 'montage')
end
% apply the following function to all the elements
% theta(x,y) = arctan(gradientX/gradientY)
divisionMatrix = arrayfun(@atan,gradientX./gradientY);
% change matrix into array
divisionArray = divisionMatrix(:);
divisionArray(isnan(divisionArray)) = [];
% convert range from (pi/2,-pi/2) to (pi,0)
for i = 1 : length(divisionArray)
     x = divisionArray(i);
     if x < 0
         divisionArray(i) = x + pi;
     end
end
% count occurence of 180 bins of angles
[counts, angles] = hist(divisionArray, 180);
% normalize all the values so that the sum is 1
counts = normalize(counts, 'norm',1);
% get kappa value for von mises distribution
% if appa  is large, the distribution becomes very concentrated to mean
% angle
kappa = circ_kappa(angles,counts);
% show the kappa value
display(kappa);
% get mean angle
[meanAngle, ul, ll] = circ_mean(divisionArray);
% fit the data into vonmises distribution(circular normal distribution)
[vonmisesP, alpha] = circ_vmpdf(angles, meanAngle, kappa);
% show polarplot
gradientsPlot = input('Would you like to create a polar plot? [y]\n','s');
if gradientsPlot == 'y';
    polarplot(alpha,vonmisesP);
    % limit the range to half the circle
    thetalim([0 180]);
end

% imported functions from CircStat2012a
function kappa = circ_kappa(alpha,w)
%
% kappa = circ_kappa(alpha,[w])
%   Computes an approximation to the ML estimate of the concentration 
%   parameter kappa of the von Mises distribution.
%
%   Input:
%     alpha   angles in radians OR alpha is length resultant
%     [w      number of incidences in case of binned angle data]
%
%   Output:
%     kappa   estimated value of kappa
%
%   References:
%     Statistical analysis of circular data, Fisher, equation p. 88
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


alpha = alpha(:);

if nargin<2
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) > size(w,1)
    w = w';
  end 
end

N = length(alpha);

if N>1
  R = circ_r(alpha,w);
else
  R = alpha;
end

if R < 0.53
  kappa = 2*R + R^3 + 5*R^5/6;
elseif R>=0.53 && R<0.85
  kappa = -.4 + 1.39*R + 0.43/(1-R);
else
  kappa = 1/(R^3 - 4*R^2 + 3*R);
end

if N<15 && N>1
  if kappa < 2
    kappa = max(kappa-2*(N*kappa)^-1,0);    
  else
    kappa = (N-1)^3*kappa/(N^3+N);
  end
end
end
function r = circ_r(alpha, w, d, dim)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
%
%   Output:
%     r		mean resultant length
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
  dim = 1;
end

if nargin < 2 || isempty(w) 
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length 
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
  c = d/2/sin(d/2);
  r = c*r;
end
end
function [mu, ul, ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit 
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
  dim = 1;
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);
% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
  t = circ_confmean(alpha,0.05,w,[],dim);
  ul = mu + t;
  ll = mu - t;
end
end
function t = circ_confmean(alpha, xi, w, d, dim)
%
% t = circ_mean(alpha, xi, w, d, dim)
%   Computes the confidence limits on the mean for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [xi   (1-xi)-confidence limits are computed, default 0.05]
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%   Output:
%     t     mean +- d yields upper/lower (1-xi)% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 5
  dim = 1;
end

if nargin < 4 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

if nargin < 3 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% set confidence limit size to default
if nargin < 2 || isempty(xi)
  xi = 0.05;
end

% compute ingredients for conf. lim.
r = circ_r(alpha,w,d,dim);
n = sum(w,dim);
R = n.*r;
c2 = chi2inv((1-xi),1);

% check for resultant vector length and select appropriate formula
t = zeros(size(r));

for i = 1:numel(r)
  if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
    t(i) = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));  % equ. 26.24
  elseif r(i) >= .9
    t(i) = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      % equ. 26.25
  else 
    t(i) = NaN;
    warning('Requirements for confidence levels not met.');
  end
end

% apply final transform
t = acos(t./R);
end
function [p, alpha] = circ_vmpdf(alpha, thetahat, kappa)

% [p alpha] = circ_vmpdf(alpha, w, p)
%   Computes the circular von Mises pdf with preferred direction thetahat 
%   and concentration kappa at each of the angles in alpha
%
%   The vmpdf is given by f(phi) =
%   (1/(2pi*I0(kappa))*exp(kappa*cos(phi-thetahat)
%
%   Input:
%     alpha     angles to evaluate pdf at, if empty alphas are chosen to
%               100 uniformly spaced points around the circle
%     [thetahat preferred direction, default is 0]
%     [kappa    concentration parameter, default is 1]
%
%   Output:
%     p         von Mises pdf evaluated at alpha
%     alpha     angles at which pdf was evaluated
%
%
%   References:
%     Statistical analysis of circular data, Fisher
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens and Marc J. Velasco, 2009
% velasco@ccs.fau.edu

% if no angles are supplied, 100 evenly spaced points around the circle are
% chosen
if nargin < 1 || isempty(alpha)
    alpha = linspace(0, 2*pi, 101)';
    alpha = alpha(1:end-1);
end
if nargin < 3
    kappa = 1;
end
if nargin < 2
    thetahat = 0;
end

alpha = alpha(:);

% evaluate pdf
C = 1/(2*pi*besseli(0,kappa));
p = C * exp(kappa*cos(alpha-thetahat));
end
         
    
         
         
    