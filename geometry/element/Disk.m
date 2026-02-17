classdef Disk < Layout2D
    %Disk: Circular planar geometry
    %
    % Purpose:
    %   Represents a disk region in a 2D plane.
    %
    % Key Properties:
    %   s - polyshape representing the disk
    %   r - The radius of the disk.
    %   o - The center coordinates [x, y] of the disk.
    %   n - The number of vertices used to approximate the circle.
    %
    % Key Methods:
    %   Disk(r, center, n)      - Construct disk geometry

    properties
        r
        o
        n
    end

    methods
        function obj = Disk(radius, center, num_points)
            arguments
                radius (1,1) {mustBePositive}
                center (1,2) double {mustBeReal} = [0,0]
                num_points (1,1) {mustBeInteger, mustBePositive} = 200
            end

            num = num_points;
            theta = linspace(0, 2*pi, num + 1); % Increase num by 1 to get num unique points
            theta(end) = []; % Remove the last point to avoid duplication at 0 and 2*pi

            coord1 = radius*cos(theta) + center(1);
            coord2 = radius*sin(theta) + center(2);

            obj@Layout2D(polyshape(coord1, coord2));
            obj.r = radius;
            obj.o = center;
            obj.n = num;
        end
    end

end