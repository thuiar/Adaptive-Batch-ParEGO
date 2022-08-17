function newxpop = MaxEI(xpop,PCheby,model,IFEs,options)
% This function is used to reccommend new multiple sample points which is with the
% best EI and optimize the MOP in the Pareto sense. 

% INPUT:
% xpop        : the observed data points
% tchebycheff : the scalar cost values/Tchebycheff function values
% model      : the GP/DACE surrogate model

% OUTPUT:
% new_xpop : new reccommending points maximizing the EI and optimizing
%               the MOP=(mu_func, sigma_func) architecture

    % Check related  parameter
    if ~isfield(options,'encoding');options.encoding = 'real';end
    if ~isfield(options,'xmin');options.xmin = zeros(1,size(xpop,1));end
    if ~isfield(options,'xmax');options.xmax = ones(1,size(xpop,1));end

    %% Maximize EI function to obtain offsprings
    offs   = [GA(xpop(TournamentSelection(2,size(xpop,1),PCheby),:));GA(xpop,{0,0,1,20})];
    N = size(offs,1); 
    EI = zeros(N,1);
    ybest = min(PCheby);
    E0 = inf;

    while IFEs>0
        for i=1:N
            [yhat,~,mse] = predictor(offs(i,:),model);
            s = sqrt(mse);

            if s==0
                EI(i)=0;
            else
               EI(i) = -(ybest-yhat)*normcdf((ybest-yhat)/s)- s*normpdf((ybest-yhat)/s);
            end
        end

        [~,index] = sort(EI);
        if EI(index(1))<E0
            best = offs(index(1),:);
            E0 = EI(index(1));
        end
        Parent = offs(index(1:ceil(N/2)),:);
        offs = [GA(Parent(TournamentSelection(2,size(Parent,1),EI(index(1:ceil(N/2)))),:));GA(Parent,{0,0,1,20})];
        IFEs = IFEs-size(offs,1);
    end
    newxpop = best;
end

