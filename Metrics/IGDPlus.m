function Score  = IGDPlus(PopObj,PF)
% <metric> <min>
% IGD+

% Written by Hongyan Wang

%------------------------------- Reference --------------------------------
% H. Ishibuchi, H. Masuda, Y. Tanigaki, and Y. Nojima. 2015. Modified distance 
% calculation in generational distance and inverted generational distance. In 
% Proceedings of the International Conference on Evolutionary Multi-Criterion 
% Optimization (EMO’15). 110–125.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    
%     if size(PopObj,2) ~= size(PF,2)
%         Score = nan;
%     else
       Distance = PF-PopObj;
       Distance(find(Distance<0))=0;
       Distance = sqrt(sum(Distance.^2,2));
       Score    = mean(Distance);
%     end
end

