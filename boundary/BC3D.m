classdef (Hidden) BC3D < BC
    %BC3D: Boundary-condition configuration for a 3D simulation domain
    %
    % Purpose:
    %   Aggregates BC settings for all three axes in a 3D model.
    %   Provides routing utilities to modify axis-specific boundary settings.
    %
    % Key Properties:
    %   bcx, bcy, bcz  - BC1D objects for x-, y-, and z-axes
    %
    % Key Methods:
    %   BC3D(bcX, bcY, bcZ)
    %       Construct from three BC1D objects.
    %
    %   setPMLThickness(axis_label, num_layer)
    %       Route to the corresponding BC1D.setPMLThickness.
    %
    %   changeType(axis_label, new_type)
    %       Route to the corresponding BC1D.changeType.

    properties (SetAccess = private)
        bcx
        bcy
        bcz
    end

    methods
        function obj = BC3D(BC_x, BC_y, BC_z)
            arguments
                BC_x BC1D
                BC_y BC1D
                BC_z BC1D
            end

            obj.bcx = BC_x;
            obj.bcy = BC_y;
            obj.bcz = BC_z;
        end

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
            switch axis_label
                case 'x'
                    targetBC = obj.bcx;
                case 'y'
                    targetBC = obj.bcy;
                case 'z'
                    targetBC = obj.bcz;
            end
        end
    end
end