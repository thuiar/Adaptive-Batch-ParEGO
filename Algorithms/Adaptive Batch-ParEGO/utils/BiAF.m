function newoffs = BiAF(xpop,tchebycheff,model,IFEs,options)
% This function is used to reccommend multiple sample points with adaptive
% selection strategy of the bi-objective acquisition function in the Pareto sense.

% The bi-objective acquisition function consists of exploitaiton (represented by mu
% function) and exploration (represented by the sqrt of sigma2 function).
% BiAF is designed to maximize the exploitation and exploration at
% the same time. It's an exploitation-exploration trade-off strategy. We
% use NSGA-II here to optimize the Sub_MOP problem in the Pareto sense.

% INPUT:
% xpop         : the observed data points
% tchebycheff : the scalar cost values/Tchebycheff function values
% model       : the GP/DACE surrogate model

% OUTPUT:
% new_offs  : new reccommending points optimizing the MOP=(mu_func, sigma_func)
%               architecture

    %% Check the related parameter
    if ~isfield(options,'ke');options.ke = 5;end
    if ~isfield(options,'alpha');options.alpha = 2;end
    if ~isfield(options,'wmax');options.wmax = 20;end
    if ~isfield(options,'NumSp');options.NumSp = 3;end
    if ~isfield(options,'epsilon');options.epsilon = 1.0;end

    [N,D] = size(xpop); 
    K = 2;
    xmin = options.xmin(1,:);
    xmax = options.xmax(1,:);
    options.ke = options.ke-1;
    CvgPar = options.epsilon;
    DivPar = 1.0 - CvgPar;

    % reference points
    [V0,options.NumSp] = UniformPoint(options.NumSp, K);
    V = V0;
    wmax = options.wmax;
    alpha = options.alpha;
    %% Generate random population
    switch options.encoding
        case 'binary'
            InitXpop = randi([0,1],N,D);
        case 'permutation'
            [~,InitXpop] = sort(rand(N,D),2); % sort the element of each row
        otherwise
            InitXpop = unifrnd(repmat(xmin,N,1),repmat(xmax,N,1));
    end

    %% Obtain objective values of initial population
    InitObjs = BiAFObj(InitXpop,xpop,tchebycheff,model,options);
    [~,~,FrontNo,CrowdDis] = EnvirSelect(InitObjs,N);

    %% Optimization Using NSGA-II
    while IFEs > 0
        MatingPool = TournamentSelection(2,N,FrontNo,-CrowdDis);
        offs  = GAEvol(InitXpop(MatingPool,:),options);
        OffsMopObjs = BiAFObj(offs,xpop,tchebycheff,model,options);
        [next,~,FrontNo,CrowdDis] = EnvirSelect([InitObjs; ...
            OffsMopObjs],N); 
        IFEs = IFEs-size(offs,1);
    end

    %% Select ke-1 nondominated solutions for GP update greedily
    MergeOffs = [InitXpop;offs];
    offs = MergeOffs(next,:);
    BiAFObjs = BiAFObj(offs,xpop,tchebycheff,model,options);
    row = size(BiAFObjs,1);

    % Scale the mu-sigma objective values to [0,1]
    NormBiAFObjs = (BiAFObjs-repmat(min(BiAFObjs,[],1),row,1)) ...
        ./repmat((max(BiAFObjs,[],1)-min(BiAFObjs,[],1)),row,1);

    %% Choose the next ke-1 solutions with adaptive selection strategy
    paras = [CvgPar,-DivPar];
    ParMatrix = repmat(paras,size(BiAFObjs,1),1);
    WeightedVars = sum(ParMatrix.*NormBiAFObjs,2); 
    [~,WeightIdx] = sort(WeightedVars);
    NewOffIdx  = WeightIdx(1:options.ke);
    newoffs = offs(NewOffIdx,:);

end