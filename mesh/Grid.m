classdef (Abstract, Hidden) Grid
    %Grid: Abstract base class for structured mesh objects (Axis/Grid2D/Grid3D)
    %
    % Key Methods (Abstract):
    %   shiftGrid(...)  - Shift grid coordinates (subclass defines signature/behavior)

    methods (Abstract)
        shiftGrid
    end

end
