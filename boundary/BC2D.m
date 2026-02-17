classdef (Hidden) BC2D < BC
    %BC2D: Boundary-condition configuration for a 2D simulation plane
    %
    % Purpose:
    %   Combines two BC1D objects for the two in-plane axes of a 2D model.
    %   Provides routing to modify per-axis BC settings through one container.
    %
    % Key Properties:
    %   bc1, bc2  - Two BC1D objects (stored in consistent cyclic axis order)
    %
    % Dependent Properties:
    %   label     - Plane label, concatenated from bc1.label and bc2.label
    %
    % Key Methods:
    %   BC2D(bcA, bcB)
    %       Construct from two BC1D objects. Requires different axis labels.
    %       Reorders inputs into cyclic axis order (x->y->z->x) for consistency.
    %
    %   setPMLThickness(axis_label, num_layer)
    %       Route to matching BC1D.setPMLThickness for the requested axis.
    %
    %   changeType(axis_label, new_type)
    %       Route to matching BC1D.changeType for the requested axis.

    properties (SetAccess = private)
        bc1
        bc2
    end

    properties (Dependent)
        label
    end

    methods
        function obj = BC2D(BC_input1, BC_input2)
            arguments
                BC_input1 BC1D
                BC_input2 BC1D
            end
            if strcmp(BC_input1.label, BC_input2.label)
                dispError('BC2D:DuplicateLabel');
            end

            % Define cyclic order: x -> y -> z -> x
            cyclic = {'x', 'y', 'z'};
            idx1 = find(strcmp(cyclic, BC_input1.label));
            idx2 = find(strcmp(cyclic, BC_input2.label));

            if mod(idx2 - idx1, 3) == 1
                obj.bc1 = BC_input1;
                obj.bc2 = BC_input2;
            else
                obj.bc1 = BC_input2;
                obj.bc2 = BC_input1;
            end
        end

        % Dependent Property
        function value = get.label(obj)
            value = [obj.bc1.label, obj.bc2.label];
        end

        % modification
        function setPMLThickness(obj, bc1d_label, num_layer)
            arguments
                obj
                bc1d_label {mustBeMember(bc1d_label, {'x','y','z'})}
                num_layer (1,1) double {mustBeInteger, mustBePositive}
            end

            % Locate and update the correct BC1D object
            targetBC = obj.findBC(bc1d_label);
            targetBC.setPMLThickness(num_layer);
        end

        function changeType(obj, bc1d_label, new_type)
            arguments
                obj
                bc1d_label {mustBeMember(bc1d_label, {'x','y','z'})}
                new_type {mustBeMember(new_type, {'Perfect', 'PML', 'Periodic'})}
            end

            % Locate and update the correct BC1D object
            targetBC = obj.findBC(bc1d_label);
            targetBC.changeType(new_type);
        end
    end

    methods (Access = private)
        function targetBC = findBC(obj, axis_label)
            if strcmp(obj.bc1.label, axis_label)
                targetBC = obj.bc1;
            elseif strcmp(obj.bc2.label, axis_label)
                targetBC = obj.bc2;
            else
                dispError('BC2D:InvalidLabel');
            end
        end
    end
end