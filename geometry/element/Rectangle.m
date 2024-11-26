classdef Rectangle < Layout2D
    % Rectangle: Vertical rectangle
    %
    % Description:
    %   The `Rectangle` class represents a rectangular shape in a 2D plane.
    %   It inherits from `Geometry2D` and provides properties and methods specific to rectangular geometries.
    %
    % Properties:
    %   u - Length unit used for the geometry (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   s - Polyshape representing the 2D shape.
    %   l1 - Length of the sides along the first dimension.
    %   l2 - Length of the sides along the second dimension.
    %
    % Methods:
    %   Rectangle(length1, length2, reference_point, ref_coord, unit) - Constructor to create an instance of the Rectangle class.
    %   printInfo() - Prints information about the rectangle object.
    %   changeUnit(new_unit) - Changes the length unit of the geometry.
    %   combineGeometry(operation, varargin) - Combines multiple Geometry2D objects based on the specified operation.
    %   dispImg(varargin) - Displays the 2D geometry.
    %   CrossSection(line) - Computes the cross-section of the geometry with a given line.
    %
    % See also:
    %   Geometry, Geometry2D

    properties
        l1
        l2
    end

    methods

         function obj = Rectangle(length1, length2, reference_point, ref_coord, unit)
            arguments
                length1 (1,1) {mustBeReal, mustBePositive}
                length2 (1,1) {mustBeReal, mustBePositive}
                reference_point {mustBeMember(reference_point,{'tl', 'bl', 'tr', 'br', 'c'})} = 'bl'
                ref_coord (1,2) {mustBeReal} = [0,0]
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            switch reference_point
                case 'tl' % top-left
                    lower_left = [ref_coord(1), ref_coord(2) - length2];
                    upper_right = [ref_coord(1) + length1, ref_coord(2)];
                case 'bl' % bottom-left
                    lower_left = [ref_coord(1), ref_coord(2)];
                    upper_right = [ref_coord(1) + length1, ref_coord(2) + length2];
                case 'tr' % top-right
                    lower_left = [ref_coord(1) - length1, ref_coord(2) - length2];
                    upper_right = [ref_coord(1), ref_coord(2)];
                case 'br' % bottom-right
                    lower_left = [ref_coord(1) - length1, ref_coord(2)];
                    upper_right = [ref_coord(1), ref_coord(2) + length2];
                case 'c' % Center
                    lower_left = [ref_coord(1)-length1/2, ref_coord(2)-length2/2];
                    upper_right = [ref_coord(1)+length1/2, ref_coord(2)+length2/2];
            end

            coord1 = [lower_left(1), upper_right(1), upper_right(1), lower_left(1)];
            coord2 = [lower_left(2), lower_left(2), upper_right(2), upper_right(2)];

            obj@Layout2D(polyshape(coord1, coord2), unit);
            obj.l1 = length1;
            obj.l2 = length2;
        end

    end

end