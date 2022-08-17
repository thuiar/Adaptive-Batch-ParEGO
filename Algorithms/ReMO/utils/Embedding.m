function A = Embedding(PDec,D,d,lower,upper,problem,options)
% This function aims to generate effective subspaces with A and Y to
% perform random embedding

% This function references ReMBO

%----------------------Reference-------------------------
% Ziyu Wang, Frank Hutter, Masrour Zoghi, David Matheson, and 
% Nando De Freitas. 2016. Bayesian optimization in a billion 
% dimensions via random embeddings. J. Artif. Int. Res. 55, 1 
% (January 2016), 361–387.
%--------------------------------------------------------

    rotate = options.rotate;
    embed = options.embed;
    ForceInBounds = options.ForceInBounds;
    if nargin < 5
        ForceInBounds = 0; 
    end

    if nargin < 5
        embed = 1;
    end

    if rotate
        % Rotate the objective function.
        [rm, ~] = qr(randn(D, D), 0);
    else
        % Do not rotate the objective function.
        rm = eye(D);
    end
    if embed
        % Generate random projection matrix A.
        [InBounds, A] = TestInBound(PDec,D,d, rm,lower,upper,problem);       
        while ~InBounds && ForceInBounds
            % Ensure that at least one maximizer fall in bound by 
            % generating as many random projection matrix A as needed.
            [InBounds, A] = TestInBound(PDec,D,d, rm,lower,upper,problem);
        end
    else
        A = eye(D, D);
    end
end

