classdef Ring < Layout2D
    %Ring: Annular (ring-shaped) planar geometry
    %
    % Purpose:
    %    Represents a ring region defined by inner and outer radii.
    %
    % Key Properties:
    %    s  - polyshape representing the annulus
    %    ri - The inner radius of the ring.
    %    ro - The outer radius of the ring.
    %    o  - The center coordinates of the ring.
    %    n  - The number of vertices used to approximate the ring.
    %
    % Key Methods:
    %    Ring([ri ro], c, n)      - Construct ring geometry
    %    combineGeometry(...)    - Combine ring with other planar geometries

    properties
        ri
        ro
        o
        n
        theta_range % Added to store the angle range
    end

    methods
        function obj = Ring(radius, center, theta_range, num_points)
            arguments
                radius (1,2) {mustBePositive}
                center (1,2) double {mustBeReal}
                theta_range (1,2) double = [0, 2*pi]
                num_points (1,1) {mustBeInteger, mustBePositive} = 200
            end

            % Ensure inner and outer radii are distinct and sorted
            if radius(1) == radius(2)
                error('Inner and outer radii must be different.');
            end
            radius = sort(radius);

            % Check if it is a full ring (360 degrees)
            isFull = abs(diff(theta_range)) >= 2*pi;

            % Generate angles based on the range provided
            % This will follow the direction of the input (e.g., [pi/2, -pi/2] goes CW)
            theta = linspace(theta_range(1), theta_range(2), num_points);

            if isFull
                % For a full ring, use the subtraction method to ensure a hole is created
                theta_full = linspace(0, 2*pi, num_points + 1);
                theta_full(end) = [];

                x_in = radius(1)*cos(theta_full) + center(1);
                y_in = radius(1)*sin(theta_full) + center(2);
                x_out = radius(2)*cos(theta_full) + center(1);
                y_out = radius(2)*sin(theta_full) + center(2);

                ring = subtract(polyshape(x_out, y_out), polyshape(x_in, y_in));
            else
                % For a segment, we create one continuous boundary
                % Trace outer arc forward, then inner arc backward
                x_out = radius(2)*cos(theta) + center(1);
                y_out = radius(2)*sin(theta) + center(2);

                x_in = radius(1)*cos(theta) + center(1);
                y_in = radius(1)*sin(theta) + center(2);

                % FLIPLR ensures the inner arc returns to the start point correctly
                ring = polyshape([x_out, fliplr(x_in)], [y_out, fliplr(y_in)]);
            end

            % Assignments
            obj@Layout2D(ring);
            obj.ri = radius(1);
            obj.ro = radius(2);
            obj.o = center;
            obj.n = num_points;
            obj.theta_range = theta_range;
        end
    end
end