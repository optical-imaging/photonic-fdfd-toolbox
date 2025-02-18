classdef Line1D < Layout1D
    % Line1D: Line in 1D space
    %
    % Description:
    %   The `Line1D` class represents a line segment in a 1D space.
    %   It inherits from `Geometry1D` and provides properties and methods specific to 1D lines.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   p1 - Coordinates of the start point of the line.
    %   p2 - Coordinates of the end point of the line.
    %
    % Methods:
    %   Line1D(point1, point2, unit) - Constructor to create an instance of the Line1D class.
    %   changeUnit(new_unit) - Changes the length unit of the geometry.
    %   getLength() - Calculates the length of the line.
    %   getCenter() - Gets the center coordinate of the line.
    %   combineGeometry(operation, varargin) - Combines multiple Line1D objects based on the specified operation.
    %
    % See also:
    %   Geometry, Geometry1D


    methods

        function obj = Line1D(point1, point2, unit)
            % Line1D constructor
            %   Constructs an instance of the Line1D class with specified start and end points and unit.
            %
            %   Syntax:
            %     obj = Line1D(point1, point2, unit)
            %
            %   Input:
            %     point1 - Coordinates of the start point (1D real values).
            %     point2 - Coordinates of the end point (1D real values).
            %     unit - Length unit (must be one of: 'm', 'cm', 'mm', 'um', 'nm', 'A'). Default is 'm'.
            %
            %   Output:
            %     obj - An instance of the Line1D class.
            
            arguments
                point1 (:,1) {mustBeReal}
                point2 (:,1) {mustBeReal}
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            obj = obj@Layout1D(point1, point2, unit);
        end

    end

end