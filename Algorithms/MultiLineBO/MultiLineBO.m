function MultiLineBO(Global)
% <algorithm> <M> <expensive>
% Multi-Line Beyesian Optimization

% IFEs --- 10000 --- Internal GA evals per iteration
    
%------------------------------- Reference --------------------------------
% [1] Kirschner J, Mutn, Mojmír, Hiller N, et al. Adaptive and Safe Bayesian
% Optimization in High Dimensions via One-Dimensional Subspaces[J]. ICML2019, 
% PMLR, Long Beach, California, USA, 2019, pp. 3429–3438.

% [2] Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------


        %% Parameter Setting
        IFEs = Global.ParameterSet(1000);
        %% Initialization
        [W,Global.N] = UniformPointNew(Global.N,Global.M);
        N             = 11*Global.D-1;
        PopDec        = UniformPointNew(N,Global.D,'Latin');
        Population    = INDIVIDUAL(repmat(Global.upper-Global.lower,N,1).*PopDec+repmat(Global.lower,N,1));
        theta         = 10;
        BestX        = UniformPointNew(1,Global.D,'Latin');
        range         = Global.upper - Global.lower;

        %% Optimization
        while Global.NotTermination(Population)
            % Preprocess
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

            % Get random direction
            direction = UniformPointNew(1,Global.D,'Latin');
            direction = direction./norm(direction);
            direction = direction .* range;

            % Get new subdomain
            temp = length(Next);
            SubDec = (Population(Next).decs- BestX(ones(1,temp),:))*(direction.');
            PCheby = PCheby(Next);

            % Eliminate the solutions having duplicated inputs or outputs
            [~,distinct1] = unique(round(SubDec*1e6)/1e6,'rows');
            [~,distinct2] = unique(round(PCheby*1e6)/1e6);
            distinct = intersect(distinct1,distinct2);
            SubDec     = SubDec(distinct,:); 
            PCheby   = PCheby(distinct);

            % Find subdomain bounds
            for i = 1:length(direction)
                if direction(i) > 0
                    TempL(i) = (Global.lower(i)-BestX(i))/direction(i);
                    TempU(i) = (Global.upper(i)-BestX(i))/direction(i);
                else
                    TempL(i) = (Global.upper(i)-BestX(i))/direction(i);
                    TempU(i) = (Global.lower(i)-BestX(i))/direction(i);
                end
            end
            Lower = max(TempL);
            Upper = min(TempU);

            % Surrogate-assisted prediction
            dmodel     = dacefit(SubDec,PCheby,'regpoly1','corrgauss',theta,1e-5,20);
            theta      = dmodel.theta;
            PopDec     = EvolEI(SubDec,PCheby,dmodel,IFEs,Lower,Upper);
            PopDec     = PopDec .* direction + BestX;
            BestX     = PopDec;
            Population = [Population,INDIVIDUAL(PopDec)];
        end
end
       