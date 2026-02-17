classdef Line2D < Layout1D
    %Line2D: 2D line-segment geometry
    %
    % Purpose:
    %   Represents line segments in a 2D plane.
    %   Used for ports and geometry cross-sections.
    %
    % Key Properties:
    %   p1  - Start point(s), size N×2
    %   p2  - End point(s),   size N×2
    %
    % Key Methods:
    %   Line2D(p1, p2)          - Construct 2D line segments

    methods
        function obj = Line2D(point1, point2)
            arguments
                point1 (:,2) {mustBeReal}
                point2 (:,2) {mustBeReal}
            end

            obj = obj@Layout1D(point1, point2);
        end
    end

end