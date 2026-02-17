classdef (Abstract, Hidden) Layout1D < Layout
    %Layout1D: Base class for 1D line-type geometry used for ports and segments
    %
    % Purpose:
    %   Stores endpoints for one or more 1D line segments (in 1D/2D/3D space).
    %   Used by Line1D/Line2D/Line3D and Port classes to define cuts/segments.
    %
    % Key Properties (protected set):
    %   p1   - Start point(s), size N×D
    %   p2   - End point(s),   size N×D
    %
    % Key Methods:
    %   Layout1D(point1, point2)      - Construct line segment(s); rejects p1==p2
    %   combineGeometry(varargin)      - (No need for 1D geometry. Put a function block to make abstract method inheritence consistant.)

    properties (SetAccess = protected)
        p1
        p2
    end

    methods
        function obj = Layout1D(point1, point2)
            arguments
                point1 double {mustBeReal}
                point2 double {mustBeReal}
            end

            % check point overlap
            if all(point1==point2)
                dispError('Geometry1D:EndpointsOverlap');
            end

            obj.p1 = point1;
            obj.p2 = point2;
        end

        function obj = combineGeometry(obj, operation, varargin)
        end
    end

end