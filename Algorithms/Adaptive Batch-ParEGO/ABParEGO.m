function ABParEGO(Global)
% <algorithm> <A><Expensive>
% Adaptive batch-ParEGO for expensive multi-objective problems
% IFEs --- 10000 --- Internal GA evals per iteration

%------------------------------- Reference --------------------------------
% [1] J. Knowles, ParEGO: A hybrid algorithm with on-line landscape
% approximation for expensive multi-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2006, 10(1): 50-66.

% [2] Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    IFEs = Global.ParameterSet(10000);
    options.ke = 5; % number of the next candidates
    options.alpha = 2; % the parameter controlling the rate of change of penalty
    options.wmax = 20; 
    options.NumSp = 4; % number of splited subspaces with weight vectors
    options.encoding = Global.encoding;
    options.xmin = Global.lower;
    options.xmax = Global.upper;
    w = 0; % inner iteration number
    InitDecay = 0.99;

    %% Generate the weight vectors and random population
	[W,Global.N] = UniformPoint(Global.N,Global.M);
    N            = 11*Global.D-1;
    PopDec       = lhsamp(N,Global.D);
    Population   = INDIVIDUAL(repmat(Global.upper-Global.lower,N,1).*PopDec+repmat(Global.lower,N,1));
	theta        = 10.*ones(1,Global.D);
    
    %% Optimization
    while Global.NotTermination(Population)
        % Decay
        options.epsilon = InitDecay^w;
        
        % Randomly select a weight vector and preprocess the data
        lamda  = W(randi(size(W,1)),:); 
        PopObj = Population.objs;
        [N,D]  = size(Population.decs);
        PopObj = (PopObj-repmat(min(PopObj,[],1),N,1))./repmat((max(PopObj,[],1)-min(PopObj,[],1)),N,1);
        PCheby = max(PopObj.*repmat(lamda,N,1),[],2)+0.05.*sum(PopObj.*repmat(lamda,N,1),2); 
        if N > 11*D-1+25
            [~,index] = sort(PCheby);
            Next      = index(1:11*D-1+25);
        else
            Next = true(N,1);
        end
        PDec   = Population(Next).decs;
        PCheby = PCheby(Next);
        
        % Eliminate the solutions having duplicated inputs or outputs
        [~,distinct1] = unique(roundn(PDec,-6),'rows');
        [~,distinct2] = unique(roundn(PCheby,-6));
        distinct = intersect(distinct1,distinct2);
        PDec     = PDec(distinct,:);
        PCheby   = PCheby(distinct);
        
        % Surrogate-assisted prediction
        dmodel     = dacefit(PDec,PCheby,'regpoly1','corrgauss',theta,1e-5.*ones(1,D),20.*ones(1,D));
        theta       = dmodel.theta;
         
        PopDecEI  = MaxEI(PDec,PCheby,dmodel,IFEs,options);
        PopDecMOP = BiAF(PDec,PCheby,dmodel,IFEs,options);
        
        % Candidate selection with Reference Vector
        % PopDecMOP = BiAFRV(PDec,PCheby,dmodel,IFEs,options);
        
        % Candidate selection with non-dominated sorting
        % PopDecMOP = BiAFNS(PDec,PCheby,dmodel,IFEs,options);
        PopDec = [PopDecEI;PopDecMOP];
        
        Population = [Population,INDIVIDUAL(PopDec)];
        w = w+1;
    end
end