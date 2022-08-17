function [neglnlike,mu_hat,sigma2_hat,U] = likelihood(x,y,theta)
% [mu,sigma2,U,neglnlike] = likelihood(x,y,theta)
% This function is used to obtain likelihood function

% INPUT:
% x          : observed data points
% y          : the Tchebycheff function values of observed data points
% theta      : correlation parameter

% OUTPUT:
% mu         : MLE (maximum likelihood estimation) mu_hat
% sigma2     : MLE sigma2_hat
% neglnlike  : MLE negtive ln likelihood function

% Make sure observed points and Tchbycheff values have the same dimensions. 
[n,m] = size(x);
ll = size(y);
if min(ll) == 1 % row vector or column vector
    y = y(:); % convert to column vector
    lenY = max(ll);
    ll = size(y);
else
    lenY = ll(1);
end
if n ~= lenY
    error('The dimension of observed points and Tchebycheff does not match!');
end

% Correlation matrix R
R = zeros(n,n);
one = ones(n,1); % one matrix
for i=1:n
    for j=i+1:n
        % upper half of correlation matrix
        R(i,j)=exp(-sum(theta.*(x(i,:)-x(j,:)).^2));
    end
end
R = R + R' + eye(n) + eye(n).*eps;

% Cholesky factorization
[U,p] = chol(R);
if p > 0 % R is not positive definite
    neglnlike = 1e4;
else % R is positive definite
%     ln_delt_R = 2*sum(log(abs(diag(U))));
    mu_hat = (one'*(U\(U'\y)))/(one'*(U\(U'\one)));
    sigma2_hat = ((y-one*mu_hat)'*(U\(U'\(y-one*mu_hat))))/n;
    neglnlike = -1*(-(n/2)*log(sigma2_hat)-0.5*log(det(R)));    
end

% max_likelihood = struct('mu',mu_hat,'sigma2',sigma2_hat,'neglnlike',neglnlike);

end

