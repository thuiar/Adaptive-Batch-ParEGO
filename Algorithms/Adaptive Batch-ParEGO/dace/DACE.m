function gpmodel = DACE(X,Y,theta0,options)
% This function is used to construct DACE model.

% INPUT:
% x          : observed points
% y          : scalar cost/Tchebycheff function values
% theta      : parameter to form correlation matrix
% options    : struct with all related paramters
% theta_low  : lower bound of theta
% theta_up   : upper bound of theta

% OUTPUT:
% gpmodel    : surrogate model, a struct

% Examine all related paprameters
if ~isfield(options, 'problem');options.problem = 'f_dtlz1';end
if ~isfield(options,'batch_size');options.batch_size = 3;end
if ~isfield(options, 'max_iter');options.max_iter = 100;end
if ~isfield(options,'min_EI');options.min_EI = 1e14;end
if ~isfield(options,'xmin');options.xmin = 0;end
if ~isfield(options,'xmax');options.xmax = 1;end
if ~isfield(options,'theta_low');options.theta_low = 1e-5.*ones(1,size(X,2));end
if ~isfield(options,'theta_up');options.theta_up = 20.*ones(1,size(X,2));end


% Make sure observed points and Tchebycheff values have the same row length 
[n,m] = size(X);
l = size(Y);
if min(l) == 1 % row vector or column vector
    Y = Y(:); % convert to column vector
    lenY = max(l);
    l = size(Y);
else
    lenY = l(1);
end
if n ~= lenY
    error('Oberved data points and Tchebycheff values must have the same row length !');
end

% Examine the upper and lower bounds of theta
theta_lower = options.theta_low;
theta_upper = options.theta_up;

len_low = size(theta_lower,2);
len_up = size(theta_upper,2);

if len_low ~=len_up || len_low ~= m ||len_up ~=m
    error('Lower bound and Upper bound of theta does not match!');
end

if any(theta_lower <=0)|| any(theta_upper<theta_lower)
    error('Lower and Upper bounds of Theta are illegal!');
end

% % Determine theta
% % Pay attention: ga is to find the minimum
% GA maximize the likelihood function
% disp('GA is beginning to maximize the likelihood function to obtain theta...');
% tic
% theta_best = ga(@(theta)likelihood(X,Y,theta),m,[],[],[],[],theta_lower,theta_upper);
% toc
% disp('GA maximizing the likelihood is completed!');

% Determine theta referenced from PlatEMO
theta_model = cal_theta(X,Y,'regpoly1','corrgauss',theta0,theta_lower,theta_upper);
theta_best = theta_model.theta;

[neglnlike,mu_hat,sigma2_hat,U]= likelihood(X,Y,theta_best);

% % return the GP surrogate model
% gpmodel = struct('mu_hat',likelihood.mu_hat,'sigma2_hat',likelihood.sigma2_hat, ...
%     'U',likelihood.U,'neglnlike',likelihood.neglnlike,'y_hat',y_hat,'s2',s2, ...
%     'mu_func',mu_func,'sigma2_func',sigma2_func,'theta',theta_best);

% return the GP surrogate model
gpmodel = struct('mu_hat',mu_hat,'sigma2_hat',sigma2_hat, ...
    'U',U,'neglnlike',neglnlike,'theta',theta_best);
end
