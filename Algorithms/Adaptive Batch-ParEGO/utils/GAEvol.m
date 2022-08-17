function offspring = GAEvol(parent,options,Parameter)
% This function is GA algorithm for evolving individuals in ParEGO

% INPUT:
% parrent    : individuals used to produce offsprings
% Parameter  : parameters used in GA

% OUTPUT:
% offspring  : evloved offsprings

% OTHERS:
% prob_corss : crossover probability
% prob_mut   : mutation probability
% dis_cross  : the distribution index of simulated binary crossover
% dis_mut    : the distribution index of polynomial mutation

%------------------------------- Reference --------------------------------
% This function is referrenced from the code of PlatEMO, "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
% https://github.com/BIMK/PlatEMO
%--------------------------------------------------------------------------

%% Parameter setting
    if nargin > 2
        [probcross,discross,probmut,dismut] = deal(Parameter{:}); % paramter distribution
    else
        [probcross,discross,probmut,dismut]=deal(1,20,1,20);
    end
    parent1 = parent(1:floor(end/2),:);
    parent2 = parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(parent1);

    encoding = options.encoding;
    switch encoding
        case 'binary'
            %% Genetic operators for binary encoding
            % One point crossover
            k = repmat(1:D,N,1) > repmat(randi(D,N,1),1,D);
            k(repmat(rand(N,1)>probcross,1,D)) = false;
            offspring1    = parent1;
            offspring2    = parent2;
            offspring1(k) = parent2(k);
            offspring2(k) = parent1(k);
            offspring     = [offspring1;offspring2];
            % Bitwise mutation
            Site = rand(2*N,D) < probmut/D;
            offspring(Site) = ~offspring(Site);
        case 'permutation'
            %% Genetic operators for permutation based encoding
            % Order crossover
            offspring = [parent1;parent2];
            k = randi(D,1,2*N);
            for i = 1 : N
                offspring(i,k(i)+1:end)   = setdiff(parent2(i,:),parent1(i,1:k(i)),'stable');
                offspring(i+N,k(i)+1:end) = setdiff(parent1(i,:),parent2(i,1:k(i)),'stable');
            end
            % Slight mutation
            k = randi(D,1,2*N);
            s = randi(D,1,2*N);
            for i = 1 : 2*N
                if s(i) < k(i)
                    offspring(i,:) = offspring(i,[1:s(i)-1,k(i),s(i):k(i)-1,k(i)+1:end]);
                elseif s(i) > k(i)
                    offspring(i,:) = offspring(i,[1:k(i)-1,k(i)+1:s(i)-1,k(i),s(i):end]);
                end
            end
        otherwise
            %% Genetic operators for real encoding
            % Simulated binary crossover
            beta = zeros(N,D);
            mu   = rand(N,D);
            beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(discross+1));
            beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(discross+1));
            beta = beta.*(-1).^randi([0,1],N,D);
            beta(rand(N,D)<0.5) = 1;
            beta(repmat(rand(N,1)>probcross,1,D)) = 1;
            offspring = [(parent1+parent2)/2+beta.*(parent1-parent2)/2
                         (parent1+parent2)/2-beta.*(parent1-parent2)/2];
            % Polynomial mutation
            Lower = repmat(options.xmin(1,:),2*N,1);
            Upper = repmat(options.xmax(1,:),2*N,1);
            Site  = rand(2*N,D) < probmut/D;
            mu    = rand(2*N,D);
            temp  = Site & mu<=0.5;
            offspring       = min(max(offspring,Lower),Upper);
            offspring(temp) = offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                              (1-(offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(dismut+1)).^(1/(dismut+1))-1);
            temp = Site & mu>0.5; 
            offspring(temp) = offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                              (1-(Upper(temp)-offspring(temp))./(Upper(temp)-Lower(temp))).^(dismut+1)).^(1/(dismut+1)));
    end
end

