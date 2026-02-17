classdef Line1D < Layout1D
    %Line1D: 1D line-segment geometry
    %
    % Purpose:
    %   Represents line segments in 1D space.
    %
    % Key Properties:
    %   p1  - Start point(s), size N×1
    %   p2  - End point(s),   size N×1
    %
    % Key Methods:
    %   Line1D(p1, p2)          - Construct 1D line segments

    methods
        function obj = Line1D(point1, point2)
            arguments
                point1 (:,1) {mustBeReal}
                point2 (:,1) {mustBeReal}
            end

            obj@Layout1D(point1, point2);
        end
    end

end