classdef Polygon < Layout2D
    %Polygon: Arbitrary polygonal planar geometry
    %
    % Purpose:
    %   Represents a closed polygon in the 2D plane.
    %
    % Key Properties:
    %   s  - polyshape representing the polygon
    %
    % Key Methods:
    %   Polygon(x, y)            - Construct polygon from vertices

    methods
        function obj = Polygon(x_coordinates, y_coordinates)
            arguments
                x_coordinates (1,:) {mustBeReal}
                y_coordinates (1,:) {mustBeReal}
            end

            obj@Layout2D(polyshape(x_coordinates, y_coordinates));
        end
    end

end