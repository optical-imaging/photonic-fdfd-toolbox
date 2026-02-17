function mustBe2D(x)
    if ~ismatrix(x)
        error('Input must be a 2D array.');
    end
end