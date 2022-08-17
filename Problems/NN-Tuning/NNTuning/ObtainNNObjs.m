function [popobj,PF]= ObtainNNObjs(X,options)
% This function is used to calculate the true objective values of X and PF of
% certain test problem

% X          : the data points
% options    : related-parameters struct
% pop_obj    : true objective function values
% PF         : true Pareto Front

if ~isfield(options,'encoding');options.encoding = 'real';end
if ~isfield(options,'NumObjs');options.NumObjs = 2;end
if ~isfield(options,'NumVars');options.NumVars = 5;end
if ~isfield(options,'popsize');options.popsize = 100;end
if ~isfield(options,'xmin');options.xmin = [1,50, -10.*ones(1,options.NumVars-2)];end
if ~isfield(options,'xmax');options.xmax = [3,500,zeros(1,options.NumVars-2)];end


problem = options.problem;
n = 10000; 
[N,M] = size(X);

for i=1:N
    [popobj(i,1),popobj(i,2)] = EvaluateModelError(X(i,:),options);
end
popobj = [popobj(:,1),popobj(:,2)];
PF = max(popobj,[],1);
end
