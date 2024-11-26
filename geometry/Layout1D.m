classdef (Abstract) Layout1D < Layout
    % Layout1D: Mid-level abstract class for 1D geometry elements
    %
    % Description:
    %   The `Layout1D` class serves as a base class for one-dimensional geometrical elements.
    %   It provides properties and methods to handle geometrical calculations, such as length,
    %   center, and unit conversion for 1D geometries. This class cannot be instantiated directly
    %   and must be subclassed.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   p1 - Coordinate of start points.
    %   p2 - Coordinate of end points.
    %
    % Methods:
    %   Layout1D(point1, point2, unit) - Constructor to create an instance of the Geometry1D class.
    %   changeUnit(unit) - Changes the length unit of the geometry.
    %   getLength() - Calculates the length of the geometry element.
    %   getCenter() - Calculates the center point of the geometry element.
    %   combineGeometry(varargin) - Combines multiple Layout1D objects.
    %
    % See also:
    %   Layout, Line1D, Line2D, Line3D

    properties (SetAccess = protected)
        p1
        p2
    end

    methods

        function obj = Layout1D(point1, point2, unit)
            % Layout1D constructor
            %   Constructs an instance of the Layout1D class with specified start and end points.
            %
            %   Syntax:
            %     obj = Layout1D(point1, point2, unit)
            %
            %   Input:
            %     point1 - Coordinates of the start points (real values).
            %     point2 - Coordinates of the end points (real values).
            %     unit - Length unit (must be one of: 'm', 'cm', 'mm', 'um', 'nm', 'A'). Default is 'm'.
            %
            %   Output:
            %     obj - An instance of the Layout1D class.

            arguments
                point1 (:,:) {mustBeReal}
                point2 (:,:) {mustBeReal}
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            % check point overlap
            if all(point1==point2)
                dispError('Geometry1D:EndpointsOverlap');
            end

            % check if points are on the same line (for Line2D and Line3D)
            if size(point1,2)>1
                Layout1D.isCollinear(point1, point2);
            end

            [p1, p2] = Layout1D.mergePoints(point1, point2);
            obj = obj@Layout(unit);
            obj.p1 = p1;
            obj.p2 = p2;
        end

        % Other methods
        function obj = changeUnit(obj, unit)
            % changeUnit - Changes the length unit of the geometry
            %
            %   Syntax:
            %     obj = changeUnit(unit)
            %
            %   Input:
            %     unit - New length unit (must be one of: 'm', 'cm', 'mm', 'um', 'nm', 'A').
            %
            %   Output:
            %     obj - An instance of the Layout1D class with updated units.

            arguments
                obj
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            obj.p1 = obj.p1*LengthUnit(unit)/LengthUnit(obj.u);
            obj.p2 = obj.p2*LengthUnit(unit)/LengthUnit(obj.u);
            obj.u = unit;
        end

        function L = getLength(obj)
            % getLength - Calculates the length of the geometry element
            %
            %   Syntax:
            %     L = getLength()
            %
            %   Output:
            %     L - Length of the geometry element.

            L = sqrt(sum((obj.p2-obj.p1).^2));
        end

        function center = getCenter(obj)
            % getCenter - Calculates the center point of the geometry element
            %
            %   Syntax:
            %     center = getCenter()
            %
            %   Output:
            %     center - Coordinates of the center point of the geometry element.

            center = (obj.p1+obj.p2)/2;
        end

        function obj = combineGeometry(obj, varargin)
            % combineGeometry - Combines multiple Layout1D objects
            %
            %   Syntax:
            %     obj = combineGeometry(geometry1, geometry2, ...)
            %
            %   Input:
            %     geometry1, geometry2, ... - Instances of Geometry1D to be combined.
            %
            %   Output:
            %     obj - Combined Layout1D object.

            arguments
                obj
            end
            arguments (Repeating)
                varargin {mustBeA(varargin, "Layout1D")}
            end

        end

    end

    methods (Static, Access = private)

        function isCollinear(p1, p2)
            % for Line2D and Line3D
            % extend 2D line to 3D expression for crossproduct
            if size(p1,2) == 2
                p1(:,3) = 0;
                p2(:,3) = 0;
            end

            % the first point is standard
            d = p2(1,:)-p1(1,:); % direction vector

            for ii = 2:size(p1,1)
                % verify parallel
                di = p2(ii,:)-p1(ii,:);
                if vecnorm(cross(d,di))>1e-10
                    dispError('Geometry1D:NotParallel', num2str(ii));
                else
                    % verify on the same line
                    qi = p1(ii,:)-p1(ii,:);
                    if vecnorm(cross(d,qi))>1e-10
                        dispError('Geometry1D:NotCollinear', num2str(ii));
                    end
                end
            end
        end

        function [p1_merge, p2_merge] = mergePoints(p1, p2)
            % initialization
            p1_merge = p1; p2_merge = p2;

            % line parametric equation
            k = p2(1,:)-p1(1,:); % [pi, pj, pk] = [bi, bj, bk] + [ki, kj, kk]*t
            nonzerok_ind = find(k~=0); % for recover later
            p1(:,k==0) = []; p2(:,k==0) = []; 
            k(k==0) = [];

            b = p1(1,:);
            t1 = (p1(:,1) - b(1))./k(1);
            t2 = (p2(:,1) - b(1))./k(1);
            [t1, t2] = Layout1D.deleteRepeatings(t1,t2);

            % chenge left and right, t2>t1
            for ii = 1:size(p1,1)
                t1i = t1(ii);
                t2i = t2(ii);
                if t2i<t1i
                    t1(ii) = t2i;
                    t2(ii) = t1i;
                end
            end

            % order points
            [t1, ind1] = sort(t1);
            t2 = t2(ind1);

            % check and merge segments
            for jj = 2:size(p1,1)
                if t1(jj)<=t2(jj-1)
                    t1(jj) = t1(jj-1);
                    t2(jj-1) = t2(jj);
                end
            end

            % delete repeat points
            flag_keep = true(size(t1));
            for ii = 2:length(t1)
                if t1(ii)==t1(ii-1) && t2(ii)==t2(ii-1)
                    flag_keep(ii) = false;
                end
            end
            [t1, t2] = Layout1D.deleteRepeatings(t1,t2);

            p1 = b+k.*t1;
            p2 = b+k.*t2;

            p1_merge = p1_merge(1:size(p1,1),:);
            p1_merge(:, nonzerok_ind) = p1;
            p2_merge = p2_merge(1:size(p2,1),:);
            p2_merge(:, nonzerok_ind) = p2;
        end

        function [t1, t2] = deleteRepeatings(t1, t2)
            flag_keep = true(size(t1));
            for ii = 2:length(t1)
                if t1(ii)==t1(ii-1) && t2(ii)==t2(ii-1)
                    flag_keep(ii) = false;
                end
            end
            t1 = t1(flag_keep);
            t2 = t2(flag_keep);
        end

    end

end