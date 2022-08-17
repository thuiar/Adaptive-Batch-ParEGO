function ReMO(Global)
% <algorithm> <R>
% Efficient global optimization for Pareto optimization. Note that ReMO is
% an optimization architecture in which the author states that ReMO can be
% equiped with any well known derivative- free MO algorithm. Since we use
% ParEGO as our baseline method, we equip ReMO with ParEGO in this paper. 

% IFEs --- 10000 --- Internal GA evals per iteration
% d --- 2 --- The number of sub-dimensions to optimize at each iteration

%------------------------------- Reference --------------------------------
% [1] H. Qian, Y. Yu, Solving high-dimensional multi-objective optimization 
% problems with low effective di- mensions, in: Proceedings of the 31st AAAI 
% Conference on Artificial Intelligence, AAAI’17, AAAI Press, 2017, p. 875–881.

% [2] Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB 
% platform for evolutionary multi-objective optimization [educational forum], 
% IEEE Computational Intelligence Magazine, 2017, 12(4): 73-87.
%--------------------------------------------------------------------------

    %% Parameter setting
    [IFEs,d] = Global.ParameterSet(10000,2);
    options.rotate = 0;
    options.ForceInBounds = 1;
    options.embed = 1;

    %% Generate the weight vectors and random population
	[W,Global.N] = UniformPoint(Global.N,Global.M);
    N            = 11*Global.D-1;
    PopDec       = lhsamp(N,Global.D);
    Population   = INDIVIDUAL(repmat(Global.upper-Global.lower,N,1).*PopDec+repmat(Global.lower,N,1));
	theta        = 10.*ones(1,Global.D);
    lower = Global.lower;
    upper = Global.upper;
    problem = Global.problem;
    
    %% Optimization
    while Global.NotTermination(Population)
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
        
        % Generate A and Y, and ensure that the Ay in the effective
        % subspace
        A =  Embedding(PDec,D,d,lower,upper,problem,options); % A=[D,d]
        PDec = PDec*A; % PDec = [n,D]*[D,d] = [n,d]
        
        % Surrogate-assisted prediction
        dmodel     = dacefit(PDec,PCheby,'regpoly1','corrgauss',theta(:,1:d),1e-5.*ones(1,d),20.*ones(1,d));
        theta      = dmodel.theta;
        TempPop = Population.decs;
        
        PopDec     = EvoALG(PCheby,TempPop*A,dmodel,IFEs,d);
        PopDec = PopDec*pinv(A); 
        Population = [Population,INDIVIDUAL(PopDec)];
    end
end