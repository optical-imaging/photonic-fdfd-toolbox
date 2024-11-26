classdef Grid2D < Grid
     % Grid2D: Even-distributed mesh plane
    %
    % Description:
    %   The `Grid2D` class represents a 2D evenly distributed mesh grid, which is a subclass of the `Grid` class.
    %   It provides methods to create, modify, and get information about 2D grids, with adjustable units and step sizes.
    %
    % Properties:
    %   axis1 - The 1st dimension axis.
    %   axis2 - The 2nd dimension axis.
    %
    % Methods:
    %   Grid2D(num_point, step_size, varargin) - Constructor to create an instance of the Grid2D class.
    %   printInfo() - Prints information about the grid properties.
    %   changeUnit(unit) - Changes the length unit of the grid.
    %   getArea(length_unit) - Calculates the area of the grid.
    %   getAxis(axis_num) - Gets one specific axis.
    %
    % See also:
    %   Grid, Axis

    properties (SetAccess = {?EigenMode2D, ?Grid3D})
        axis1
        axis2
    end

    methods

        % Constructor
        function obj = Grid2D(num_point,step_size,varargin)
            % Grid2D constructor
            %   Constructs an instance of the Grid2D class with specified number of points, step sizes, and additional options.
            %
            %   Syntax:
            %     obj = Grid2D(num_point, step_size, 'Center', center, 'Unit', unit, 'ZeroPoint', zero_point)
            %
            %   Input:
            %     num_point - Number of points along each axis (2-element vector of positive integers).
            %     step_size - Step size between points along each axis (2-element vector of positive real values).
            %     Name-Value arguments:
            %       'Center' - Center value of the grid for each axis (default is [0, 0]).
            %       'Unit' - Length unit for the grid (default is 'm').
            %       'ZeroPoint' - Indicator if each axis includes zero (default is [1, 1]).
            %
            %   Output:
            %     obj - An instance of the Grid2D class.

            arguments
                num_point (1,2) double {mustBePositive, mustBeInteger}
                step_size (1,2) double {mustBeReal,mustBePositive}
            end

            arguments(Repeating)
                varargin
            end

            % Set up inputParser to handle name-value pairs
            p = inputParser;
            addParameter(p, 'Center', [0, 0], @(x) isnumeric(x) && numel(x) == 2); % Default center is 0
            addParameter(p, 'Unit', 'm', @(x) ismember(x,...
                {'m', 'cm', 'mm', 'um', 'nm', 'A'})); % Default unit is 'm'
            addParameter(p, 'ZeroPoint', [1, 1], @(x) isnumeric(x) && all(ismember(x, [0, 1])) && numel(x) == 2);
            parse(p, varargin{:});
            center = p.Results.Center;
            length_unit = p.Results.Unit;
            zero_point = p.Results.ZeroPoint;
            
            obj.axis1 = Axis(num_point(1),step_size(1), 'Center', center(1), 'Unit', length_unit, 'ZeroPoint', zero_point(1));
            obj.axis2 = Axis(num_point(2),step_size(2), 'Center', center(2), 'Unit', length_unit, 'ZeroPoint', zero_point(2));
        end

        % Display
        function printInfo(obj)
            info_name = {'Object name';
                         'Number of points';
                         'Length unit'
                         'Step size';
                         'Range';
                         'Total area size'};
            value = {inputname(1); 
                     [num2str(obj.axis1.n),' x ',num2str(obj.axis2.n)];
                     obj.u;
                     [num2str(obj.axis1.d),' x ',num2str(obj.axis2.d)];
                     ['[',num2str(obj.axis1.v(1)),' , ',num2str(obj.axis1.v(end)),...
                        '] x [',num2str(obj.axis2.v(1)),' , ',num2str(obj.axis2.v(end)),']'];
                     [num2str(obj.axis1.getLength),' x ',num2str(obj.axis2.getLength)]};
            
            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength + 1, info_name{ii}, value{ii});
            end
        end

        % Other methods
        function obj = changeUnit(obj,unit)
            % changeUnit - Changes the length unit of the grid
            %
            %   Syntax:
            %     obj = changeUnit(unit)
            %
            %   Input:
            %     unit - The new length unit for the grid.
            %
            %   Output:
            %     obj - The Grid2D object with updated length unit for both axes.

            arguments
                obj
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})}
            end
            
            obj.axis1 = obj.axis1.changeUnit(unit);
            obj.axis2 = obj.axis2.changeUnit(unit);
        end
        
        function [S,area_unit] = getArea(obj,length_unit)
            % getArea - Calculates the area of the grid
            %
            %   Syntax:
            %     [S, area_unit] = getArea(length_unit)
            %
            %   Input:
            %     length_unit - Desired length unit for the output (default is 'm').
            %
            %   Output:
            %     S - Total area of the grid in the specified unit.
            %     area_unit - The unit used for the area output.

            arguments
                obj 
                length_unit {mustBeMember(length_unit,{'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            L1 = obj.axis1.getLength(length_unit);
            L2 = obj.axis2.getLength(length_unit);
            S = L1*L2;
            area_unit = [length_unit,'^2'];
        end

        function axis = getAxis(obj, axis_num)
            % getAxis - Gets one specific axis of the grid
            %
            %   Syntax:
            %     axis = getAxis(axis_num)
            %
            %   Input:
            %     axis_num - The axis number to retrieve (1 or 2).
            %
            %   Output:
            %     axis - The specified axis of the grid.
            
            arguments
                obj
                axis_num {mustBeMember(axis_num, [1,2])}
            end

            switch axis_num
                case 1
                    axis = obj.axis1;
                case 2
                    axis = obj.axis2;
            end
        end

    end

end