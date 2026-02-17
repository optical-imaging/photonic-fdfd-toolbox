classdef Grid2D < Grid
    %Grid2D: 2D uniform grid defined by two orthogonal Axis objects
    %
    % Key Properties:
    %   axis1                  - First axis (Axis object; label matches plane_label(1))
    %   axis2                  - Second axis (Axis object; label matches plane_label(2))
    %   label                  - Plane label 'xy'|'yz'|'zx' (dependent; [axis1.label axis2.label])
    %
    % Key Methods:
    %   Grid2D(N2, d2, label)        - Construct 2D grid from [Nx Ny], [dx dy], and plane label
    %   get.label()                  - Return plane label from axis labels
    %   shiftGrid(dim, mode, val)    - Shift one axis ('x'|'y'|'z') by min/max/center alignment
    %   getAxis(axis_id)             - Return axis by index (1|2) or label ('x'|'y'|'z')

    properties (SetAccess = {?EigenMode2D, ?Grid3D})
        axis1
        axis2
    end

    properties (Dependent)
        label
    end

    methods
        % constructor
        function obj = Grid2D(num_point, step_size, plane_label)
            arguments
                num_point (1,2) double {mustBePositive, mustBeInteger}
                step_size (1,2) double {mustBeReal, mustBePositive}
                plane_label string {mustBeMember(plane_label, {'xy','yz','zx'})}
            end

            plane_label = char(plane_label); % "xy" is also allowed 

            obj.axis1 = Axis(num_point(1), step_size(1), plane_label(1));
            obj.axis2 = Axis(num_point(2), step_size(2), plane_label(2));
        end

        % dependent
        function val = get.label(obj)
            val = [obj.axis1.label, obj.axis2.label];
        end

        % modify values of axis
        function obj = shiftGrid(obj, dim, mode, new_value)
            arguments
                obj
                dim {mustBeMember(dim, {'x','y','z'})}
                mode {mustBeMember(mode, {'min','max','c'})}
                new_value (1,1) double {mustBeReal}
            end

            grid_label = {obj.axis1.label, obj.axis2.label};
            match_idx = find(strcmp(dim, grid_label), 1);

            if isempty(match_idx)
                dispError('Grid2D:AxisNotIncluded', obj.label, obj.label(1), obj.label(2));
            end

            switch match_idx
                case 1
                    obj.axis1 = obj.axis1.shiftGrid(mode, new_value);
                case 2
                    obj.axis2 = obj.axis2.shiftGrid(mode, new_value);
            end
        end

        % extract one axis
        function axis = getAxis(obj, axis_id)
            arguments
                obj
                axis_id {mustBeMember(axis_id, {1, 2, 'x', 'y', 'z'})}
            end
            
            % Case 1: Numeric index (1 or 2)
            if isnumeric(axis_id)
                switch axis_id
                    case 1
                        axis = obj.axis1;
                    case 2
                        axis = obj.axis2;
                end
            else
                % Case 2: Label input ('x', 'y', or 'z')
                label_id = char(axis_id);
                
                if strcmp(label_id, obj.axis1.label)
                    axis = obj.axis1;
                elseif strcmp(label_id, obj.axis2.label)
                    axis = obj.axis2;
                else
                    dispError('Grid2D:AxisNotIncluded', obj.label, obj.label(1), obj.label(2));
                end
            end
        end
    end
end
