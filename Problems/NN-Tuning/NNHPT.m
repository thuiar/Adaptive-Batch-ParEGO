classdef NNHPT < PROBLEM
% <problem> <NNHPT>
% Hyper-parameter tunning of neural network

%----------------------------------Note-----------------------------------
% The hyper-parameter tuning task is referenced to the sourcecode of study [1]. 
% We thank the author very much for his contribution, and the Github repository
% of the source code is: https://github.com/rasmusbergpalm/DeepLearnToolbox

%------------------------------- Reference --------------------------------
% [1]Palm R B. Prediction as a candidate for learning deep hierarchical models of 
% data. Technical University of Denmark, 2012, 5.

% [2] E. Zitzler, K. Deb, and L. Thiele, Comparison of multiobjective
% evolutionary algorithms: Empirical results, Evolutionary computation,
% 2000, 8(2): 173-195.

% [3] Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform 
% for evolutionary multi-objective optimization [educational forum], IEEE  
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%  X range:
%  xmin = [1 50  0   -10 -10 -10 3 800 0 0.4]
%  xmax = [3 500 0.9 0  0  0 6 1200 2 0.6]

    methods
        %% Initialization
        function obj = NNHPT()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = [1,50,0, -10.*ones(1,obj.Global.D-8),1, 800, 0, 0.4, 0];
            obj.Global.upper    = [3,500,0.9,zeros(1,obj.Global.D-8),3, 1500, 2, 0.6,1];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            row = size(PopDec,1);
            for i=1:row
                plot = 0;
                if i<4
                    plot = 0;
                end
                [PopObj(i,1),PopObj(i,2)] = EvaluateModelError(PopDec(i,:),plot);
            end
            PopObj = [PopObj(:,1),PopObj(:,2)];
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = [1 10];
        end
    end
end

