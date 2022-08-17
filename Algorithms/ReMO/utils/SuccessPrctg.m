function [prct,A] = SuccessPrctg(D, num_trial, dim,d, bounds, ParetoSet, rm)
    total = 0;
    cmbnts = nchoosek(1:d,dim); % matrix-all possible combinations of a set of values
    num_PS = size(ParetoSet, 2);

    for i = 1:num_trial
        indices = 1:D;
        A = randn(D, d);

        if nargin > 6
            A = rm*A;
        end
        fail = 1;

        for j = 1:size(cmbnts, 1)
            for k = 1:num_PS
                YFound = inv(A(indices(1:dim), cmbnts(j, :))) *  ParetoSet(:, k); % A-1x* = y*
                if ~(sum(YFound <= bounds(:,2)) < dim || sum(YFound >= bounds(:,1)) < dim)
                    fail = 0;
                end
            end
        end
        total = total + fail;
    end
    prct = 1 - total/num_trial; % 0/1
end

