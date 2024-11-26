classdef Grid3D < Grid
    % Grid3D: Even-distributed mesh volume
    %
    % Description:
    %   The `Grid3D` class represents a 3D evenly distributed mesh grid, which is a subclass of the `Grid` class.
    %   It provides methods to create, modify, and get information about 3D grids, with adjustable units and step sizes.
    %
    % Properties:
    %   axisx - The x-axis.
    %   axisy - The y-axis.
    %   axisz - The z-axis.
    %
    % Methods:
    %   Grid3D(num_point, step_size, varargin) - Constructor to create an instance of the Grid3D class.
    %   printInfo() - Prints information about the grid properties.
    %   changeUnit(unit) - Changes the length unit of the grid.
    %   getVolume(length_unit) - Calculates the volume of the grid.
    %   getAxis(axis_label) - Gets one specific axis of the grid.
    %   getPlane(plane_label) - Gets one specific plane of the grid.
    %
    % See also:
    %   Grid, Axis, Grid2D

    properties (SetAccess = private)
        axisx
        axisy
        axisz
    end

    methods

        function obj = Grid3D(num_point,step_size,varargin)
            % Grid3D constructor
            %   Constructs an instance of the Grid3D class with specified number of points, step sizes, and additional options.
            %
            %   Syntax:
            %     obj = Grid3D(num_point, step_size, 'Center', center, 'Unit', unit, 'ZeroPoint', zero_point)
            %
            %   Input:
            %     num_point - Number of points along each axis (3-element vector of positive integers).
            %     step_size - Step size between points along each axis (3-element vector of positive real values).
            %     Name-Value arguments:
            %       'Center' - Center value of the grid for each axis (default is [0, 0, 0]).
            %       'Unit' - Length unit for the grid (default is 'm').
            %       'ZeroPoint' - Indicator if each axis includes zero (default is [1, 1, 1]).
            %
            %   Output:
            %     obj - An instance of the Grid3D class.

            arguments
                num_point (1,3) double {mustBePositive, mustBeInteger}
                step_size (1,3) double {mustBeReal,mustBePositive}
            end

            arguments(Repeating)
                varargin
            end

            % Set up inputParser to handle name-value pairs
            p = inputParser;
            addParameter(p, 'Center', [0, 0, 0], @(x) isnumeric(x) && numel(x) == 3); % Default center is 0
            addParameter(p, 'Unit', 'm', @(x) ismember(x,...
                {'m', 'cm', 'mm', 'um', 'nm', 'A'})); % Default unit is 'm'
            addParameter(p, 'ZeroPoint', [1, 1, 1], @(x) isnumeric(x) && all(ismember(x, [0, 1])) && numel(x) == 3);

            parse(p, varargin{:});
            center = p.Results.Center;
            length_unit = p.Results.Unit;
            zero_point = p.Results.ZeroPoint;
            
            obj.axisx = Axis(num_point(1),step_size(1), 'Center', center(1), 'Unit', length_unit, 'ZeroPoint', zero_point(1));
            obj.axisy = Axis(num_point(2),step_size(2), 'Center', center(2), 'Unit', length_unit, 'ZeroPoint', zero_point(2));
            obj.axisz = Axis(num_point(3),step_size(3), 'Center', center(3), 'Unit', length_unit, 'ZeroPoint', zero_point(3));
        end

        % Display
        function printInfo(obj)
            info_name = {'Object name';
                'Number of points';
                'Length unit';
                'Step size';
                'Range';
                'Total volume size'};
            value = {inputname(1);
                [num2str(obj.axisx.n),' x ',num2str(obj.axisy.n),' x ',num2str(obj.axisz.n)];
                obj.u;
                [num2str(obj.axisx.d),' x ',num2str(obj.axisy.d),' x ',num2str(obj.axisz.d)];
                ['[',num2str(obj.axisx.v(1)),' ,',num2str(obj.axisx.v(end)),...
                '] x [',num2str(obj.axisy.v(1)),' ,',num2str(obj.axisy.v(end)),']',...
                '] x [',num2str(obj.axisz.v(1)),' ,',num2str(obj.axisz.v(end)),']'];
                [num2str(obj.axisx.getLength),' x ',num2str(obj.axisy.getLength),' x '...
                ,num2str(obj.axisz.getLength)]};

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength + 1, info_name{ii}, value{ii});
            end
        end

        % Other methods
        function [V,unit] = getVolume(obj,length_unit)
            % getVolume - Calculates the volume of the grid
            %
            %   Syntax:
            %     [V, unit] = getVolume(length_unit)
            %
            %   Input:
            %     length_unit - Desired length unit for the output (default is 'm').
            %
            %   Output:
            %     V - Total volume of the grid in the specified unit.
            %     unit - The unit used for the volume output.

            arguments
                obj 
                length_unit {mustBeMember(length_unit,{'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            Lx = obj.axisx.getLength(length_unit);
            Ly = obj.axisy.getLength(length_unit);
            Lz = obj.axisz.getLength(length_unit);
            V = Lx*Ly*Lz;
            unit = [length_unit,'^3'];
        end

        function axis = getAxis(obj, axis_label)
            % getAxis - Gets one specific axis of the grid
            %
            %   Syntax:
            %     axis = getAxis(axis_label)
            %
            %   Input:
            %     axis_label - The axis label to retrieve ('x', 'y', or 'z').
            %
            %   Output:
            %     axis - The specified axis of the grid.

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

        function plane = getPlane(obj, plane_label)
            % getPlane - Gets one specific plane of the grid
            %
            %   Syntax:
            %     plane = getPlane(plane_label)
            %
            %   Input:
            %     plane_label - The plane label to retrieve ('xy', 'yz', or 'zx').
            %
            %   Output:
            %     plane - The specified plane of the grid as a Grid2D object.

            arguments
                obj
                plane_label char {mustBeMember(plane_label, {'xy','yz','zx'})}
            end

            plane = Grid2D(1,1);
            switch plane_label
                case 'xy'
                    plane.axis1 = obj.axisx;
                    plane.axis2 = obj.axisy;
                case 'yz'
                    plane.axis1 = obj.axisy;
                    plane.axis2 = obj.axisz;
                case 'zx'
                    plane.axis1 = obj.axisz;
                    plane.axis2 = obj.axisx;
            end
        end

    end

end