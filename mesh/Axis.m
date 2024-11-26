classdef Axis < Grid
    % Axis: Coordinate axis
    %
    % Description:
    %   The `Axis` class represents an even-distributed coordinate axis, which is a subclass of the `Grid` class.
    %   It provides methods to create, modify, and get information about coordinate axes, with adjustable units and step sizes.
    %
    % Properties:
    %   v - Coordinate axis data.
    %   u - Length unit for the axis (e.g., 'm', 'cm', 'mm', 'um', 'nm', 'A').
    %   n - Number of points along the axis (dependent property).
    %   d - Step size between points along the axis (dependent property).
    %
    % Methods:
    %   Axis(num_point, step_size, varargin) - Constructor to create an instance of the Axis class.
    %   printInfo() - Prints information about the axis properties.
    %   getLength(length_unit) - Calculates the total length of the axis.
    %   changeUnit(new_unit) - Changes the length unit of the axis.
    %
    % See also:
    %   Grid

    properties (SetAccess = {?Port1D, ?Port2D})
        v 
        u
    end

    properties (Dependent)
        n
        d
    end

    methods

        % Constructor
        function obj = Axis(num_point,step_size,varargin)
             % Axis constructor
            %   Constructs an instance of the Axis class with specified number of points, step size, and additional options.
            %
            %   Syntax:
            %     obj = Axis(num_point, step_size, 'Center', center, 'Unit', unit, 'ZeroPoint', zero_point)
            %
            %   Input:
            %     num_point - Number of points along the axis (positive integer).
            %     step_size - Step size between points along the axis (positive real value).
            %     Name-Value arguments:
            %       'Center' - Center value of the axis (default is 0).
            %       'Unit' - Length unit for the axis (default is 'm').
            %       'ZeroPoint' - Indicator if axis includes zero (default is 1).
            %
            %   Output:
            %     obj - An instance of the Axis class.

            arguments
                num_point (1,1) double {mustBePositive, mustBeInteger}
                step_size (1,1) double {mustBeReal,mustBePositive}
            end

            arguments(Repeating)
                varargin
            end

            % Set up inputParser to handle name-value pairs
            p = inputParser;
            addParameter(p, 'Center', 0, @(x) isnumeric(x) && isscalar(x)); % Default center is 0
            addParameter(p, 'Unit', 'm', @(x) ismember(x,...
                {'m', 'cm', 'mm', 'um', 'nm', 'A'})); % Default unit is 'm'
            addParameter(p, 'ZeroPoint', 1, @(x) ismember(x,[0,1]));
            parse(p, varargin{:});
            center = p.Results.Center;
            length_unit = p.Results.Unit;
            zero_point = p.Results.ZeroPoint;

            % Even number of points
            if mod(num_point,2) == 0 && zero_point == 1
                obj.v = center + ((-floor((num_point-1)/2)):(ceil((num_point-1)/2)))*step_size;
            else
                obj.v = center + (-(num_point-1)/2:(num_point-1)/2)*step_size;
            end

            % Get value
            obj.u = length_unit;
        end

        % Dependent
        function value = get.n(obj)
            value = numel(obj.v);
        end

        function value = get.d(obj)
            value = abs(obj.v(1)-obj.v(2));
        end

        % Display
        function printInfo(obj)
            info_name = {'Axis name';
                'Number of points';
                'Length unit'
                'Step size';
                'Axis range';
                'Total Length'};
            value = {inputname(1);
                num2str(obj.n);
                obj.u;
                num2str(obj.d);
                ['[',num2str(obj.v(1)),',',num2str(obj.v(end)),']'];
                num2str(obj.getLength)};

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength + 1, info_name{ii}, value{ii});
            end
        end

        % Other methods
        function [L,unit] = getLength(obj,length_unit)
            % getLength - Calculates the total length of the axis
            %
            %   Syntax:
            %     [L, unit] = getLength(length_unit)
            %
            %   Input:
            %     length_unit - Desired length unit for the output (default is the current unit of the axis).
            %
            %   Output:
            %     L - Total length of the axis in the specified unit.
            %     unit - The unit used for the length output.

            arguments
                obj
                length_unit {mustBeMember(length_unit,{'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            if nargin == 1
                unit = obj.u;
            else
                unit = length_unit;
            end

            L = obj.d*LengthUnit(obj.u)*obj.n/LengthUnit(unit);
        end

        function obj = changeUnit(obj, new_unit)
            % changeUnit - Changes the length unit of the axis
            %
            %   Syntax:
            %     obj = changeUnit(new_unit)
            %
            %   Input:
            %     new_unit - The new length unit for the axis.
            %
            %   Output:
            %     obj - The Axis object with updated length unit and coordinate data.
            
            arguments
                obj
                new_unit {mustBeMember(new_unit,{'m', 'cm', 'mm', 'um', 'nm', 'A'})}
            end

            obj.v = obj.v * LengthUnit(obj.u) / LengthUnit(new_unit);
            obj.u = new_unit;  % Update the unit property
        end

    end

    methods (Access = {?Device, ?FDFD, ?Source})

        function obj = cutAxis(obj, v_idx)
            obj.v = obj.v(v_idx);
        end

        function obj = doublesampleAxis(obj)
            N_ds = 2*obj.n;
            d_ds = 1/2*obj.d;
            center_ds = (obj.v(1)+obj.v(end))/2-d_ds;

            obj = Axis(N_ds, d_ds, "Unit", obj.u, "Center", center_ds, "ZeroPoint", 0);
        end
        
    end

end