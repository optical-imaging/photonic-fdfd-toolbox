classdef Grid3D < Grid
    %Grid3D: 3D uniform grid volume defined by three orthogonal Axis objects
    %
    % Key Properties:
    %   axisx                     - x-axis (Axis)
    %   axisy                     - y-axis (Axis)
    %   axisz                     - z-axis (Axis)
    %
    % Key Methods:
    %   Grid3D(N3, d3)                 - Construct 3D grid from [Nx Ny Nz], [dx dy dz]
    %   shiftGrid(dim, mode, val)      - Shift one axis ('x'|'y'|'z') by min/max/center alignment
    %   getAxis(axis_label)            - Return axis by label 'x'|'y'|'z'
    %   getPlane(plane_label)          - Return Grid2D plane: 'xy'|'yz'|'zx' (axes reused, not copied)

    properties (SetAccess = private)
        axisx
        axisy
        axisz
    end

    methods
        % constructor
        function obj = Grid3D(num_point, step_size)
            arguments
                num_point (1,3) double {mustBePositive, mustBeInteger}
                step_size (1,3) double {mustBeReal, mustBePositive}
            end

            obj.axisx = Axis(num_point(1), step_size(1), 'x');
            obj.axisy = Axis(num_point(2), step_size(2), 'y');
            obj.axisz = Axis(num_point(3), step_size(3), 'z');
        end

        % modify grid value
        function obj = shiftGrid(obj, dim, mode, new_value)
            arguments
                obj
                dim {mustBeMember(dim, {'x','y','z'})}
                mode {mustBeMember(mode, {'min','max','c'})}
                new_value (1,1) double {mustBeReal}
            end
            switch dim
                case 'x'
                    obj.axisx = obj.axisx.shiftGrid(mode, new_value);
                case 'y'
                    obj.axisy = obj.axisy.shiftGrid(mode, new_value);
                case 'z'
                    obj.axisz = obj.axisz.shiftGrid(mode, new_value);
            end
        end

        % extract one axis
        function axis = getAxis(obj, axis_label)
            arguments
                obj
                axis_label char {mustBeMember(axis_label, {'x','y','z'})}
            end
            switch axis_label
                case 'x'
                    axis = obj.axisx;
                case 'y'
                    axis = obj.axisy;
                case 'z'
                    axis = obj.axisz;
            end
        end

        % extract one plane
        function plane = getPlane(obj, plane_label)
            arguments
                obj
                plane_label string {mustBeMember(plane_label, {'xy', 'yz', 'zx'})}
            end
            
            switch plane_label
                case 'xy'
                    plane = Grid2D([obj.axisx.n, obj.axisy.n], [obj.axisx.d, obj.axisy.d], 'xy');
                    plane.axis1 = obj.axisx; 
                    plane.axis2 = obj.axisy;
                case 'yz'
                    plane = Grid2D([obj.axisy.n, obj.axisz.n], [obj.axisy.d, obj.axisz.d], 'yz');
                    plane.axis1 = obj.axisy;
                    plane.axis2 = obj.axisz;
                case 'zx'
                    plane = Grid2D([obj.axisz.n, obj.axisx.n], [obj.axisz.d, obj.axisx.d], 'zx');
                    plane.axis1 = obj.axisz;
                    plane.axis2 = obj.axisx;
            end
        end
    end
end
