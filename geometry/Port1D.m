classdef Port1D < Line2D
    %Port1D: 1D port defined as a line segment in a 2D simulation plane
    %
    % Purpose:
    %   Represents a source/monitor port for 2D and Semi-3D simulations.
    %
    % Key Properties:
    %   p1     - Port start point, size 1×2
    %   p2     - Port end point,   size 1×2
    %   label  - Plane label associated with the port
    %   dir    - Injection direction ('+x','-x','+y','-y')
    %
    % Key Methods:
    %   Port1D(...)             - Construct axis-aligned port definition

    properties
        label % working plane
        dir
    end

    methods
        function obj = Port1D(plane_label, direction, position, span_range)
            arguments
                plane_label {mustBeMember(plane_label, {'xy', 'yz', 'zx'})}
                direction {mustBeMember(direction, {'+x', '-x', '+y', '-y', '+z', '-z'})}
                position {mustBeReal}
                span_range (1,2) {mustBeReal}
            end

            % check label
            if ~contains(plane_label, direction(2))
                dispError('Port1D:PortLableNotMatch', direction(2), plane_label);
            end

            % check span range
            if span_range(1)==span_range(2)
                dispError('Port1D:RangeSidesEqual');
            end
            span_range = sort(span_range);

            % construct port endpoints
            ind_axis_port = find(char(plane_label)==direction(2));
            switch ind_axis_port
                case 1
                    point1 = [position, span_range(1)];
                    point2 = [position, span_range(2)];
                case 2
                    point1 = [span_range(1), position];
                    point2 = [span_range(2), position];
            end

            obj@Line2D(point1, point2);
            obj.label = plane_label;
            obj.dir = direction;
        end
    end

end