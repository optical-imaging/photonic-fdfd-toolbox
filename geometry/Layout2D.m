classdef (Abstract) Layout2D < Layout
    % Layout2D: Mid-level abstract class for 2D geometry elements
    %
    % Description:
    %   The `Layout2D` class serves as a base class for two-dimensional geometrical elements.
    %   It provides properties and methods to handle geometrical operations, such as combining shapes,
    %   changing units, and visualization. This class cannot be instantiated directly and must be subclassed.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   s - Polyshape object representing the 2D geometry.
    %
    % Methods:
    %   Layout2D(poly_shape, unit) - Constructor to create an instance of the Layout2D class.
    %   combineGeometry(operation, varargin) - Combines multiple Layout2D objects based on the specified operation.
    %   dispImg(varargin) - Displays the 2D geometry.
    %   changeUnit(new_unit) - Changes the length unit of the geometry.
    %   CrossSection(line) - Computes the cross-section of the geometry with a given line.
    %
    % See also:
    %   Layout, Disk, Rectangle, Polygon, Ring

    properties (SetAccess = protected)
        s
    end

    methods

        function obj = Layout2D(poly_shape, unit)
            arguments
                poly_shape polyshape
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            obj@Layout(unit);
            obj.s = poly_shape;
        end

        function obj = combineGeometry(obj, operation, varargin)
            arguments
                obj
                operation {mustBeMember(operation, {'Union', 'Intersect', 'Subtract', 'ExclusiveOR'})}
            end
            arguments (Repeating)
                varargin
            end

            % check and unify units
            [varargin{:}] = obj.unifyUnit(varargin{:});

            switch operation
                case 'Union'
                    for ii = 1:numel(varargin)
                        obj.s = union(obj.s, varargin{ii}.s);
                    end
                case 'Intersect'
                    for ii = 1:numel(varargin)
                        obj.s = intersect(obj.s, varargin{ii}.s);
                    end
                case 'Subtract'
                    for ii = 1:numel(varargin)
                        obj.s = subtract(obj.s, varargin{ii}.s);
                    end
                case 'ExclusiveOR'
                    for ii = 1:numel(varargin)
                        obj.s = xor(obj.s, varargin{ii}.s);
                    end
            end
        end

        % Display
        function dispImg(obj, varargin)
            figure;
            plot(obj.s, varargin{:});
            axis image;
        end

        % Other methods
        function obj = changeUnit(obj, new_unit)
            arguments
                obj
                new_unit {mustBeMember(new_unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})}
            end

            vertice = (original_shape.Vertices)*LengthUnit(new_unit)/LengthUnit(obj.u);
            obj.s = polyshape(vertice(:, 1), vertice(:, 2));
            obj.u = new_unit;
        end

        function cross_section = CrossSection(obj, line)
            arguments
                obj
                line {mustBeA(line, "Line2D")}
            end

            line.changeUnit(obj.u);

            [bi, bj] = boundary(obj.s);
            segment_intersect1 = []; segment_intersect2 = [];
            for ii = size(line.p1,1)
                % find intersection points
                [intersect_pointi, intersect_pointj] =...
                    polyxpoly([line.p1(ii,1), line.p2(ii,1)], [line.p1(ii,2), line.p2(ii,2)], bi, bj);
                key_point = [line.p1(ii,:); [intersect_pointi, intersect_pointj]; line.p2(ii,:)];
                t = vecnorm(key_point - line.p1(ii,:), 2, 2);
                [~, idx] = sort(t);
                key_point = key_point(idx, :);
                % find midpoints inside polyshpae
                for jj = 1:size(key_point,1)-1
                    midpoint = (key_point(jj,:)+key_point(jj+1,:))/2;
                    if isinterior(obj.s, midpoint(1), midpoint(2))
                        segment_intersect1 = [segment_intersect1; key_point(jj,:)];
                        segment_intersect2 = [segment_intersect2; key_point(jj+1,:)];
                    end
                end

                cross_section = Line2D(segment_intersect1, segment_intersect2, obj.u);
            end
        end

    end

    methods (Access = private, Static)

        function varargout = unifyUnit(varargin)
            arguments (Repeating)
                varargin {mustBeA(varargin, "Layout2D")}
            end

            % check
            units = cellfun(@(x) x.u, varargin, 'UniformOutput', false);
            if all(strcmp(units, units{1}))
                varargout = varargin;
            else
                % find majority unit
                [uniqueUnits, ~, unitIndices] = unique(units);
                counts = accumarray(unitIndices, 1);
                [~, majorityIndex] = max(counts);
                majorityUnit = uniqueUnits{majorityIndex};

                dispWarning('Geometry2D:UnitsNotAlign', majorityUnit);

                for ii = 1:numel(varargin)
                    varargout{ii} = varargin{ii}.changeUnit(majorityUnit);
                end
            end
        end

    end
    
end