classdef Axis < Grid
    %Axis: 1D uniform coordinate axis used as a mesh primitive and 1D grid
    %
    % Key Properties:
    %   v              - Coordinate vector (centered by construction)
    %   label          - Axis label: 'x' | 'y' | 'z'
    %   n              - Number of grid points (dependent; numel(v))
    %   d              - Grid spacing (dependent; |v(2)-v(1)|)
    %
    % Key Methods:
    %   Axis(n, d, label)       - Construct a centered uniform axis with n points and spacing d
    %   get.n()                 - Return number of points
    %   get.d()                 - Return spacing
    %   shiftGrid(mode, val)    - Shift axis so min/max/center aligns to target value
    %   cutAxis(idx)            - Keep only v(idx) (internal trimming used by Device/FDFD/Source)
    %   doublesampleAxis()      - Return 2×-sampled axis (internal for 2× maps / Yee sampling)

    properties (SetAccess = {?Port1D, ?Port2D})
        v
        label = 'x'
    end

    properties (Dependent)
        n
        d
    end

    methods
        % constructor
        function obj = Axis(num_point, step_size, axis_label)
            arguments
                num_point (1,1) double {mustBePositive, mustBeInteger}
                step_size (1,1) double {mustBeReal, mustBePositive}
                axis_label char {mustBeMember(axis_label, {'x', 'y', 'z'})}
            end

            obj.v = ((-floor((num_point-1)/2)):(ceil((num_point-1)/2))) * step_size;
            % Note: For odd number of points, the central point is 0 as
            % default; for even number of points, the central-left point is
            % 0 as default

            obj.label = char(axis_label);
        end

        % dependent properties
        function value = get.n(obj)
            value = numel(obj.v);
        end

        function value = get.d(obj)
            value = abs(obj.v(1) - obj.v(2));
        end

        % modify values of axis
        function obj = shiftGrid(obj, mode, new_value)
            arguments
                obj
                mode {mustBeMember(mode, {'min','max','c'})}
                new_value (1,1) double {mustBeReal}
            end

            switch mode
                case 'min'
                    shift_val = new_value - obj.v(1);
                case 'max'
                    shift_val = new_value - obj.v(end);
                case 'c'
                    shift_val = new_value - (obj.v(1)+obj.v(end))/2;
            end

            obj.v = obj.v + shift_val;
        end
    end

    methods (Access = {?Device, ?FDFD, ?Source})

        function obj = cutAxis(obj, v_idx)
            obj.v = obj.v(v_idx);
        end

        function axis_ds = doublesampleAxis(obj)
            N_ds = 2 * obj.n; % is an even number
            d_ds = 0.5 * obj.d;

            axis_ds = Axis(N_ds, d_ds, obj.label); % even number points axis
            axis_ds = axis_ds.shiftGrid('min', obj.v(1));
        end
    end
end