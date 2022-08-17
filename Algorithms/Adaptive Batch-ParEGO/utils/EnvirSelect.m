function [Next,Population,FrontNo,CrowdDis] = EnvirSelect(Population,N)
% The environmental selection of NSGA-II
% This code is referenced from the PlatEMO
% N : the number of individuals to be sorted at least. 
% The individuals have not been sorted are assigned a front number of inf.

%------------------------------- Reference --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true; %200
    
    %% Population for next generation
    Population = Population(Next,:);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end