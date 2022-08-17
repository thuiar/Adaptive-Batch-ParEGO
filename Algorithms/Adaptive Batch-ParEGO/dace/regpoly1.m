function  [f, df] = regpoly1(S)
%REGPOLY1  First order polynomial regression function
%
% Call:    f = regpoly1(S)
%          [f, df] = regpoly1(S)
%
% S : m*n matrix with design sites
% f = [1  s]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update April 12, 2002

[m n] = size(S);
f = [ones(m,1)  S]; % 在S前拼接一列单位元素
if  nargout > 1
  df = [zeros(n,1) eye(n)]; % 在n维单位矩阵前拼接一列0元素
end