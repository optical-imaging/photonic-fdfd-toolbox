classdef (Hidden) BC1D < BC
    %BC1D: Boundary-condition configuration for a single coordinate axis
    %
    % Purpose:
    %   Stores boundary type and PML thickness for one axis ('x','y','z').
    %   Used by solvers to decide boundary enforcement on that axis.
    %
    % Key Properties (read-only):
    %   label   - Axis label: 'x' | 'y' | 'z'
    %   type    - Boundary type: 'Perfect' | 'PML' | 'Periodic'
    %   pmlnum  - PML thickness in grid cells (0 unless type == 'PML')
    %
    % Key Methods:
    %   BC1D(axis_label, BC_type)
    %       Construct BC for a given axis and type. If type=='PML', pmlnum
    %       initializes to default (currently 8), else pmlnum=0.
    %
    %   setPMLThickness(num_layers)
    %       Set pmlnum (only allowed when type == 'PML').
    %
    %   changeType(new_type)
    %       Update type. If new_type=='PML', pmlnum resets to default (8);
    %       otherwise pmlnum resets to 0.

    properties (SetAccess = private)
        label
        type
        pmlnum
    end

    methods
        function obj = BC1D(axis_label, BC_type)
            arguments
                axis_label {mustBeMember(axis_label, {'x','y','z'})}
                BC_type {mustBeMember(BC_type, {'Perfect', 'PML', 'Periodic'})}
            end
            obj.label = axis_label;
            obj.type = BC_type;

            if strcmp(obj.type, 'PML')
                obj.pmlnum = 8; % Default to 8 if not specified
            else
                obj.pmlnum = 0; % Force pmlnum to 0 for Perfect and Periodic
            end
        end

        function setPMLThickness(obj, num_layers)
            arguments
                obj
                num_layers (1,1) double {mustBeInteger, mustBePositive}
            end
            if ~strcmp(obj.type, 'PML')
                dispError('BC1D:PMLSetOnly');
            end
            obj.pmlnum = num_layers;
        end

        function changeType(obj, new_type)
            arguments
                obj
                new_type {mustBeMember(new_type, {'Perfect', 'PML', 'Periodic'})}
            end
            obj.type = new_type;

            % refresh numbers of layers of PML
            if strcmp(obj.type, "PML")
                obj.pmlnum = 8;
            else
                obj.pmlnum = 0;
            end
        end
    end
end