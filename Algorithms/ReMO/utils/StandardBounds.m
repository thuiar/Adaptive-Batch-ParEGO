function bounds = StandardBounds(dim)
    bounds = ones(dim, 2); 
    bounds(:, 1) = -bounds(:, 1);
end