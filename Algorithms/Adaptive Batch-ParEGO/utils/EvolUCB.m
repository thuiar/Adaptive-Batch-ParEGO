function new_xpop = EvolUCB(xpop,tchebycheff,model,options)
% This function is used to reccommend new multiple sample points which is with the
% best EI and optimize the MOP in the Pareto sense. 

% INPUT:
% xpop        : the observed data points
% tchebycheff : the scalar cost values/Tchebycheff function values
% model       : the GP/DACE surrogate model

% OUTPUT:
% new_xpop    : new reccommending points maximizing the EI and optimizing
%               the MOP=(mu_func, sigma_func) architecture

    % Check related  parameter
    if ~isfield(options,'encoding');options.encoding = 'real';end
    if ~isfield(options,'max_eval');options.max_eval = 10000;end
    if ~isfield(options,'xmin');options.xmin = zeros(1,size(xpop,1));end
    if ~isfield(options,'xmax');options.xmax = ones(1,size(xpop,1));end

    max_eval = options.max_eval;
    %% Maximize EI function to obtain offsprings
    % Tournament selection for xpop to select size(xpop,1) individuals, they
    % may be reproductive
    offs = [GA_EVOL(xpop(TournamentSelection(2,size(xpop,1),tchebycheff),:),options); ...
        GA_EVOL(xpop,options,{0,0,1,20})];
    N = size(offs,1); % number of offsprings
    LCB = zeros(N,1); % pre-allocation
    ybest = min(tchebycheff);
    LCB0 = inf;

    while max_eval>0
        for i=1:N
            [yhat,~,mse] = predictor(offs(i,:),model);
            s = sqrt(mse);
            LCB(i) = yhat - 2*s;
        end

        [~,index] = sort(LCB);
        if LCB(index(1))<LCB0
            best = offs(index(1),:);
            LCB0 = LCB(index(1));
        end
        parent = offs(index(1:ceil(N/2)),:);
        offs = [GA_EVOL(parent(TournamentSelection(2,size(parent,1),LCB(index(1:ceil(N/2)))),:), ...
            options); GA_EVOL(parent,options,{0,0,1,20})];
        max_eval = max_eval-size(offs,1);
    end
    new_xpop = best;
end

