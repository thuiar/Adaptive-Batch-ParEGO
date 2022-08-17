function [InBound,A] = TestInBound(PDec,D,d,rm,lower,upper,problem)
% Check whether Ay is in the effective subspace

% This function references ReMBO

%----------------------Reference-------------------------
% Ziyu Wang, Frank Hutter, Masrour Zoghi, David Matheson, and 
% Nando De Freitas. 2016. Bayesian optimization in a billion 
% dimensions via random embeddings. J. Artif. Int. Res. 55, 1 
% (January 2016), 361–387.
%--------------------------------------------------------

    ParetoSet = GetPS(PDec,d,lower,upper,problem);
    TestBounds = StandardBounds(2);
    scale = max(1.5*log(d));
    TestBounds = TestBounds*scale;

    [prct, A] = SuccessPrctg(D, 1, 2, d, TestBounds,...
        ParetoSet, rm);
    InBound =  prct;
end

