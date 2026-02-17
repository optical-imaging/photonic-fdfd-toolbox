classdef Line3D < Layout1D
    %Line3D: 3D line-segment geometry
    %
    % Purpose:
    %   Represents line segments in 3D space.
    %
    % Key Properties:
    %   p1  - Start point(s), size N×3
    %   p2  - End point(s),   size N×3
    %
    % Key Methods:
    %   Line3D(p1, p2)          - Construct 3D line segments

    methods
        function obj = Line3D(point1, point2)
            arguments
                point1 (:,3) {mustBeReal}
                point2 (:,3) {mustBeReal}
            end

            obj = obj@Layout1D(point1, point2);
        end
    end

end