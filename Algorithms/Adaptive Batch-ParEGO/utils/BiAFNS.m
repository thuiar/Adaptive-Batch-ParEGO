function NewOffs = BiAFNS(xpop,tchebycheff,model,IFEs,options)
% This function is used to reccommend multiple sample points with non-dominated
% sorting strategy of the bi-objective acquisition function in the Pareto sense. 

% The bi-objective acquisition function consists of exploitaiton (represented by mu
% function) and exploration (represented by the sqrt of sigma2 function).
% It is designed to maximize the exploitation and exploration at
% the same time, which is an exploitation-exploration trade-off strategy. We
% use NSGA-II here to optimize the Sub_MOP problem in the Pareto sense.

% INPUT:
% xpop         : the observed data points
% tchebycheff : the scalar cost values/Tchebycheff function values
% model       : the GP/DACE surrogate model

% OUTPUT:
% new_offs  : new reccommending points optimizing the MOP=(mu_func, sigma_func)
%               architecture

    %% Check the related parameter
    if ~isfield(options,'encoding');options.encoding = 'real';end
    if ~isfield(options,'popsize');options.popsize = 100;end
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
    N = options.popsize; 
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
            [~,InitXpop] = sort(rand(N,D),2);
        otherwise
            InitXpop = unifrnd(repmat(xmin,N,1),repmat(xmax,N,1));
    end
    %% Obtain objective values of init_xpops
    InitMopObjs = BiAFObj(InitXpop,xpop,tchebycheff,model,options);
    [~,~,FrontNo,CrowdDis] = EnvirSelect(InitMopObjs,N);

    %% Optimization Using NSGAII
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

    %% k-means cluster
    rng(1); % For reproducibility
    cindex = kmeans(BiAFObjs,20);
    OffIdx = [];

    for i=1:20 % choose one solution from each cluster
        index = find(cindex == i);
        ChosenIdx = index(randi(length(index)));
        OffIdx = [OffIdx,ChosenIdx];   
    end

    NewOffs = offs(OffIdx,:);
    NewBiAFObjs = BiAFObjs(OffIdx,:);

    NonObjs = NDSort(NewBiAFObjs,N)==1;
    NondomiOffs = offs(NonObjs==1,:);
    NonBiAFObjs = NewBiAFObjs(NonObjs==1,:);

    if size(NondomiOffs,1) == options.ke
%         disp ('The first Pareto front is equal to ke!');
        NewOffs = NondomiOffs;
    elseif size(NondomiOffs,1) > options.ke
%         disp ('The first Pareto front is enough!');
        NewOffs = RandChoose(NondomiOffs,NonBiAFObjs,options); 
    elseif size(NondomiOffs,1) < options.ke
%         disp ('The first Pareto front is not enough!');
        NewOffs = ChooseMore(offs,NondomiOffs,NewBiAFObjs,N,options);
    end

end

%% k-means clustering to select ke solutions from the first Pareto front
function NewOffs = RandChoose(offs,NonObjs,options)
    OffIdx = randperm(size(offs,1),options.ke);
    NewOffs = offs(OffIdx,:);
end

%% choose more solution from higher Pareto front
function NewOffs = ChooseMore(offs,NondomiOffs,objs,N,options)

res = options.ke - size(NondomiOffs,1);
[FrontNum,MaxFront]= NDSort(objs,N);
ResIdx = [];

for i=2:MaxFront
    if res>0
        domi = NDSort(objs,N)==i;
        DomiOffs = offs(domi==1,:);
        if size(DomiOffs,1)> res
            DomiIdx = find(domi);
            ResIdx = [ResIdx,DomiIdx(randperm(size(DomiOffs,1),res))];
        elseif size(DomiOffs,1) < res
            ResIdx = [ResIdx,find(domi)];
        elseif size(DomiOffs,1) == res
            ResIdx = [ResIdx,find(domi)];
        end 
        lengthres = length(ResIdx);
        res = res -length(ResIdx);
    end
end

MoreOffs = offs(ResIdx,:);
NewOffs = [NondomiOffs;MoreOffs];
end

