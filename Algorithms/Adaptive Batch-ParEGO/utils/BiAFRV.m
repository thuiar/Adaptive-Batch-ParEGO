function NewOffs = BiAFRV(xpop,tchebycheff,model,IFEs,options)
% This function is used to reccommend multiple sample points with reference vector
% of the bi-objective acquisition function in the Pareto sense. 
 
% The bi-objective acquisition function consists of exploitaiton (represented by mu
% function) and exploration (represented by the sqrt of sigma2 function).
% It is designed to maximize the exploitation and exploration at
% the same time, which is an exploitation-exploration trade-off strategy. We
% use NSGA-II here to optimize the Sub_MOP problem in the Pareto sense.

% INPUT:
% xpop        : the observed data points
% tchebycheff : the scalar cost values/Tchebycheff function values
% model       : the GP/DACE surrogate model

% OUTPUT:
% new_offs        : new reccommending points optimizing the MOP=(mu_func, sigma_func)
%               architecture

%% Check the related parameter
if ~isfield(options,'encoding');options.encoding = 'real';end
if ~isfield(options,'PopSize');options.PopSize = 100;end
if ~isfield(options,'xmin');options.xmin = zeros(1,size(xpop,1));end
if ~isfield(options,'xmax');options.xmax = ones(1,size(xpop,1));end
if ~isfield(options,'ke');options.ke = 5;end
if ~isfield(options,'alpha');options.alpha = 2;end
if ~isfield(options,'wmax');options.wmax = 20;end
if ~isfield(options,'NumSp');options.NumSp = 3;end
%--------------------------------------------------------------------------
% N is the number of individuals to be sorted at least. N=1 indicates
% finding only the first non-dominated front, N = size(xpop,1)/2 indicates
% sorting only half the population (which is often used in the algorithm).
N = options.PopSize; 
K = 2; 
D = size(xpop,2);
xmin = options.xmin(1,:);
xmax = options.xmax(1,:);
options.ke = options.ke-1;

% reference points
[~,options.NumSp] = UniformPoint(options.NumSp, K);
%% Generate random population
switch options.encoding
    case 'binary'
        InitXpop = randi([0,1],N,D);
    case 'permutation'
        [~,InitXpop] = sort(rand(N,D),2); % sort the element of each row
    otherwise
        % Generate 100 initial population
        InitXpop = unifrnd(repmat(xmin,N,1),repmat(xmax,N,1));
end

%% Obtain objective values of InitXpops
InitMopObjs = BiAFObj(InitXpop,xpop,tchebycheff,model,options);
[~,~,FrontNo,CrowdDis] = EnvirSelect(InitMopObjs,N);

%% Optimization Using NSGA-II
while IFEs > 0
    MatingPool = TournamentSelection(2,N,FrontNo,-CrowdDis);
    offs  = GAEvol(InitXpop(MatingPool,:),options);
    OffsMopObjs = BiAFObj(offs,xpop,tchebycheff,model,options);
    [next,~,FrontNo,CrowdDis] = EnvirSelect([InitMopObjs; ...
        OffsMopObjs],N); 
    IFEs = IFEs-size(offs,1);
end

%% Select ke-1 nondominated solutions for GP update using Reference vector
MergeOffs = [InitXpop;offs];
offs = MergeOffs(next,:);
BiAFObjs = BiAFObj(offs,xpop,tchebycheff,model,options);
objs = BiAFObjs;

IdxRF = CandiSelectRV(objs,options);
if length(IdxRF)>= options.ke
    offs = offs(IdxRF,:);
    NewOffIdx = randperm(size(offs,1),options.ke);
else
    NewOffIdx = [IdxRF,randperm(size(offs,1),options.ke-length(IdxRF))];
end

NewOffs = offs(NewOffIdx,:);

end