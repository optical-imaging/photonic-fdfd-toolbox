classdef Line3D < Layout1D
    % Line3D: Line in 3D space
    %
    % Description:
    %   The `Line3D` class represents a line segment in a 3D space.
    %   It inherits from `Geometry1D` and provides properties and methods specific to 3D lines.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   p1 - Coordinates of the start point of the line.
    %   p2 - Coordinates of the end point of the line.
    %
    % Methods:
    %   Line3D(point1, point2, unit) - Constructor to create an instance of the Line3D class.
    %   changeUnit(new_unit) - Changes the length unit of the geometry.
    %   getLength() - Calculates the length of the line.
    %   getCenter() - Gets the center coordinate of the line.
    %   combineGeometry(operation, varargin) - Combines multiple Line1D objects based on the specified operation.
    %
    % See also:
    %   Geometry, Geometry1D

    methods

        function obj = Line3D(point1, point2, unit)
            arguments
                point1 (:,3) {mustBeReal}
                point2 (:,3) {mustBeReal}
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            obj = obj@Layout1D(point1, point2, unit);
        end

    end

end