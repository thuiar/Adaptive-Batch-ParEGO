function index = CandiSelectRV(PopObj,options)
% This function is used to locate the ke-1 candidate points in the
% bi-objective acqusition function space for GP updation in a greedy mode

if ~isfield(options,'ke');options.ke = 5;end
if ~isfield(options,'alpha');options.alpha = 2;end
if ~isfield(options,'wmax');options.wmax = 20;end
if ~isfield(options,'NumSp');options.NumSp = 3;end

%% Initialize the refence vector to split the space into two sub-spaces
[N,K] = size(PopObj);
[V0,options.NumSp] = UniformPoint(options.NumSp, K);
V = V0(2,:)/norm(V0(2)); % Unitization
V(:,2) = -V(:,2);

%% Translate the population
PopObj = PopObj - repmat(min(PopObj,[],1),N,1);


%% Locate candidates with better predictive values in sub-spaces
DotPopV = dot(PopObj',repmat(V,N,1)');
index = find(DotPopV <= 0);

end