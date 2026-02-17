function mustBeRealAndGE1(x)
    if ~isreal(x)
        error('All elements must be real.');
    end
    if any(x(:) < 1)
        error('All elements must be greater than or equal to 1.');
    end
end