classdef Ring < Layout2D
    % Ring: Concentric ring
    %
    % Description:
    %   The `Ring` class represents a concentric ring in a 2D plane.
    %   It inherits from `Geometry2D` and provides properties and methods specific to ring geometries.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   s - Polyshape representing the 2D shape.
    %   ri - Inner radius of the ring.
    %   ro - Outer radius of the ring.
    %   o - Center coordinates of the ring.
    %   n - Number of vertices used to approximate the ring.
    %
    % Methods:
    %   Ring(radius, center, unit, num_points) - Constructor to create an instance of the Ring class.
    %   printInfo() - Prints information about the ring object.
    %   changeUnit(new_unit) - Changes the length unit of the geometry.
    %   combineGeometry(operation, varargin) - Combines multiple Geometry2D objects based on the specified operation.
    %   dispImg(varargin) - Displays the 2D geometry.
    %   CrossSection(line) - Computes the cross-section of the geometry with a given line.
    %
    % See also:
    %   Geometry, Geometry2D

    properties
        ri
        ro
        o
        n
    end

    methods

        function obj = Ring(radius, center, unit, num_points)
            % Ring constructor
            %   Constructs an instance of the Ring class with specified inner and outer radii, center, unit, and number of points.
            %
            %   Syntax:
            %     obj = Ring(radius, center, unit, num_points)
            %
            %   Input:
            %     radius - A two-element vector specifying the inner and outer radii (positive values).
            %     center - Center coordinates of the ring (real values).
            %     unit - Length unit (must be one of: 'm', 'cm', 'mm', 'um', 'nm', 'A'). Default is 'm'.
            %     num_points - Number of vertices used to approximate the ring (positive integer). Default is 200.
            %
            %   Output:
            %     obj - An instance of the Ring class.
            
            arguments
                radius (1,2) {mustBePositive}
                center (1,2) double {mustBeReal}
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
                num_points (1,1) {mustBeInteger, mustBePositive} = 200
            end

            % check inner and outer radius
            if radius(1)==radius(2)
                dispError('Ring:RingRadiusEqual');
            else
                radius = sort(radius);
            end

            num = num_points;
            theta = linspace(0, 2*pi, num+1);
            theta(end) = [];

            x_in = radius(1)*cos(theta) + center(1);
            y_in = radius(1)*sin(theta) + center(2);
            disk_in = polyshape(x_in, y_in);

            x_out = radius(2)*cos(theta) + center(1);
            y_out = radius(2)*sin(theta) + center(2);
            disk_out = polyshape(x_out, y_out);

            ring = subtract(disk_out, disk_in);

            obj@Layout2D(ring, unit);
            obj.ri = radius(1);
            obj.ro = radius(2);
            obj.o = center;
            obj.n = num_points;
        end

        % Display
        function printInfo(obj)
            info_name = {'Object Name';
                'Unit';
                'Inner radius';
                'Outer radius';
                'Center'};
            value = {inputname(1);
                obj.u;
                num2str(obj.ri);
                num2str(obj.ro);
                num2str(['[',num2str(obj.c(1)),',',num2str(obj.c(2)),']'])};

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength+1, info_name{ii}, value{ii});
            end
        end

    end

end