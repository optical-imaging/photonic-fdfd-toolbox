classdef Polygon < Layout2D
    % Polygon: Closed polygon
    %
    % Description:
    %   The `Polygon` class represents a closed polygon in a 2D plane.
    %   It inherits from `Geometry2D` and provides properties and methods specific to polygon geometries.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   s - Polyshape representing the 2D shape.
    %
    % Methods:
    %   Polygon(x_coordinates, y_coordinates, unit) - Constructor to create an instance of the Polygon class.
    %   printInfo() - Prints information about the polygon object.
    %   changeUnit(new_unit) - Changes the length unit of the geometry.
    %   combineGeometry(operation, varargin) - Combines multiple Geometry2D objects based on the specified operation.
    %   dispImg(varargin) - Displays the 2D geometry.
    %   CrossSection(line) - Computes the cross-section of the geometry with a given line.
    %
    % See also:
    %   Geometry, Geometry2D

methods

    function obj = Polygon(x_coordinates, y_coordinates, unit)
            % Polygon constructor
            % Polygon constructor
            %   Constructs an instance of the Polygon class with specified x and y coordinates and unit.
            %
            %   Syntax:
            %     obj = Polygon(x_coordinates, y_coordinates, unit)
            %
            %   Input:
            %     x_coordinates - X-coordinates of the polygon vertices (real values).
            %     y_coordinates - Y-coordinates of the polygon vertices (real values).
            %     unit - Length unit (must be one of: 'm', 'cm', 'mm', 'um', 'nm', 'A'). Default is 'm'.
            %
            %   Output:
            %     obj - An instance of the Polygon class.
            
        arguments
            x_coordinates (1,:) {mustBeReal}
            y_coordinates (1,:) {mustBeReal}
            unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
        end

        obj@Layout2D(polyshape(x_coordinates, y_coordinates), unit);
    end
    
end

end