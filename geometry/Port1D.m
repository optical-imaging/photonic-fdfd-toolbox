classdef Port1D < Line2D
    % Port1D: 1D port defined in 2D space
    %
    % Description:
    %   The `Port1D` class represents a 1D port defined in a 2D space.
    %   It inherits from `Line2D` and provides properties and methods specific to ports.
    %   The port can only be placed either vertically or horizontally in the plane.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   p1 - Coordinates of the start point of the port.
    %   p2 - Coordinates of the end point of the port.
    %   label - Plane label for identifying the port's plane.
    %   dir - Direction of the port.
    %
    % Methods:
    %   Port1D(plane_label, direction, position, span_range, unit) - Constructor to create an instance of the Port1D class.
    %   changeUnit(new_unit) - Changes the length unit of the geometry.
    %   getLength() - Calculates the length of the port.
    %   getCenter() - Gets the center coordinate of the port.
    %   combineGeometry(operation, varargin) - Combines multiple Line1D objects based on the specified operation.
    %
    % See also:
    %   Geometry, Geometry1D, Line2D

    properties
        dir
    end

    methods
        function obj = Port1D(plane_label, direction, position, span_range, unit)
            % Port1D constructor
            %   Constructs an instance of the Port1D class with specified plane label, direction, position, span range, and unit.
            %
            %   Syntax:
            %     obj = Port1D(plane_label, direction, position, span_range, unit)
            %
            %   Input:
            %     plane_label - Plane label for identifying the port's plane ('xy', 'yz', 'zx').
            %     direction - Direction of the port ('+x', '-x', '+y', '-y', '+z', '-z').
            %     position - Position value of the port (real value).
            %     span_range - Span range along the corresponding plane axis (two-element vector).
            %     unit - Length unit (must be one of: 'm', 'cm', 'mm', 'um', 'nm', 'A'). Default is 'm'.
            %
            %   Output:
            %     obj - An instance of the Port1D class.
            
            arguments
                plane_label {mustBeMember(plane_label, {'xy', 'yz', 'zx'})}
                direction {mustBeMember(direction, {'+x', '-x', '+y', '-y', '+z', '-z'})}
                position {mustBeReal}
                span_range (1,2) {mustBeReal}
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            % check span range
            if span_range(1)==span_range(2)
                dispError('Port1D:RangeSidesEqual');
            end
            span_range = sort(span_range);

            % construct port endpoints
            ind_axis_port = find(plane_label==direction(2));
            switch ind_axis_port
                case 1
                    point1 = [position, span_range(1)];
                    point2 = [position, span_range(2)];
                case 2
                    point1 = [span_range(1), position];
                    point2 = [span_range(2), position];
            end

            obj@Line2D(point1, point2, unit);
            obj.dir = direction;
        end

    end

end