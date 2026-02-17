classdef Rectangle < Layout2D
    %Rectangle: Axis-aligned rectangular planar geometry
    %
    % Purpose:
    %   Represents a rectangle in the 2D plane.
    %
    % Key Properties:
    %   s  - polyshape representing the rectangle
    %   l1 - The side length along the first dimension (e.g., x).
    %   l2 - The side length along the second dimension (e.g., y).
    %
    % Key Methods:
    %   Rectangle(l1, l2, ref, coord) - Construct rectangle using reference convention

    properties
        l1
        l2
    end

    methods
        function obj = Rectangle(length1, length2, reference_point, ref_coord)
            arguments
                length1 (1,1) {mustBeReal, mustBePositive}
                length2 (1,1) {mustBeReal, mustBePositive}
                reference_point {mustBeMember(reference_point,{'tl', 'bl', 'tr', 'br', 'c'})} = 'c'
                ref_coord (1,2) {mustBeReal} = [0,0]
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

            obj@Layout2D(polyshape(coord1, coord2));
            obj.l1 = length1;
            obj.l2 = length2;
        end
    end

end