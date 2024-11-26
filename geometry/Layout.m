classdef (Abstract) Layout
    % Layout: Top-level abstract class for the Geometry family
    %
    % Description:
    %   The `Layout` class serves as a base class for different types of geometries.
    %   It provides a common interface for handling geometrical objects and operations.
    %   This abstract class cannot be instantiated directly and must be subclassed.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %
    % Abstract Methods:
    %   changeUnit - Abstract method that should be implemented to change the unit of the geometry.
    %   combineGeometry - Abstract method that should be implemented to combine geometries.
    %
    % Methods:
    %   Layout(unit) - Constructor to create an instance of the Geometry class.
    %
    % See also:
    %   Layout1D, Layout2D


    properties (SetAccess = protected)
        u
    end

    methods (Abstract)
        changeUnit
        combineGeometry
    end

    methods

        function obj = Layout(unit)
            % Layout constructor
            %   Constructs an instance of the Geometry class with a specified length unit.
            %
            %   Syntax:
            %     obj = Layout(unit)
            %
            %   Input:
            %     unit - A string specifying the length unit. Must be one of the following:
            %            'm', 'cm', 'mm', 'um', 'nm', 'A'. Default is 'm'.
            %
            %   Output:
            %     obj - An instance of the Geometry class (to be subclassed).

            arguments
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end
            obj.u = unit;
        end

    end

end