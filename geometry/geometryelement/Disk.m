classdef Disk < Layout2D
    % Disk: Circle disk
    %
    % Description:
    %   The `Disk` class represents a circular disk in a 2D plane.
    %   It inherits from `Geometry2D` and provides properties and methods specific to circular geometries.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   s - Polyshape representing the 2D shape.
    %   r - Radius of the disk.
    %   o - Center coordinates of the disk.
    %   n - Number of vertices used to approximate the disk.
    %
    % Methods:
    %   Disk(radius, center, unit, num_points) - Constructor to create an instance of the Disk class.
    %   printInfo() - Prints information about the disk object.
    %   changeUnit(new_unit) - Changes the length unit of the geometry.
    %   combineGeometry(operation, varargin) - Combines multiple Geometry2D objects based on the specified operation.
    %   dispImg(varargin) - Displays the 2D geometry.
    %   CrossSection(line) - Computes the cross-section of the geometry with a given line.
    %
    % See also:
    %   Geometry, Geometry2D

    properties
        r
        o
        n
    end

    methods

        function obj = Disk(radius, center, unit, num_points)
            % Disk constructor
            %   Constructs an instance of the Disk class with specified radius, center, length unit, and number of points.
            %
            %   Syntax:
            %     obj = Disk(radius, center, unit, num_points)
            %
            %   Input:
            %     radius - Radius of the disk (positive value).
            %     center - Center coordinates of the disk (real values).
            %     unit - Length unit (must be one of: 'm', 'cm', 'mm', 'um', 'nm', 'A'). Default is 'm'.
            %     num_points - Number of vertices used to approximate the disk (positive integer). Default is 200.
            %
            %   Output:
            %     obj - An instance of the Disk class.
            
            arguments
                radius (1,1) {mustBePositive}
                center (1,2) double {mustBeReal}
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
                num_points (1,1) {mustBeInteger, mustBePositive} = 200
            end

            num = num_points;
            theta = linspace(0, 2*pi, num + 1); % Increase num by 1 to get num unique points
            theta(end) = []; % Remove the last point to avoid duplication at 0 and 2*pi
            
            coord1 = radius*cos(theta) + center(1);
            coord2 = radius*sin(theta) + center(2);
            
            obj@Layout2D(polyshape(coord1, coord2), unit);
            obj.r = radius;
            obj.o = center;
            obj.n = num;
        end

        % Display
        function printInfo(obj)
            info_name = {'Object Name';
                'Unit';
                'Radius';
                'Center'};
            value = {inputname(1);
                obj.u;
                num2str(obj.r);
                num2str(['[',num2str(obj.c(1)),',',num2str(obj.c(2)),']'])};

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength+1, info_name{ii}, value{ii});
            end
        end

    end

end