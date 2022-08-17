function objs =BiAFObj(InitXpop,xpop,tchebycheff,model,options)
% This function is used to calculate the objective values of bi-objective 
% acquisition function of given population xpop

% INPUT:
% options    : struct of related optimizing parameters

% OUTPUT:
% objs       : Sub-MOP objevtive values

    % Check related parameter
    if ~isfield(options,'encoding');options.encoding = 'real';end
    if ~isfield(options,'xmin');options.xmin = zeros(1,size(xpop,1));end
    if ~isfield(options,'xmax');options.xmax = ones(1,size(xpop,1));end

    objs = zeros(size(InitXpop,1),2); % Two sub-objectives in bi-objective acquisition function
    [MuFunc,Sigma2Func] = predictor(InitXpop,model);
    SigmaFunc = sqrt(Sigma2Func);
    objs=[MuFunc,-SigmaFunc];

end

