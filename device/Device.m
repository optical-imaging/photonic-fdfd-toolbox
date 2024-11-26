classdef (Abstract) Device < handle
    % Device: Abstract base class for devices in the simulation.
    %
    % Description:
    %   The `Device` class provides a framework for defining different devices in a simulation.
    %   It handles device geometry, material properties, meshing, and unit changes.
    %   This abstract class cannot be instantiated directly and must be subclassed.
    %   The `handle` inheritance allows instances of the `Device` class to be passed by reference,
    %   which means that modifications to an object are reflected across all references to that object.
    %
    % Properties:
    %   geo - Geometry of the device.
    %   mat - Material properties of the device.
    %   layer - Layer number of the device.
    %   meshflag - Flag indicating if the device has been meshed.
    %   mesh - Meshing data for the device.
    %   eps - Relative permittivity of the device.
    %   mu - Relative permeability of the device.
    %
    % Abstract Methods:
    %   meshDevice() - Abstract method that should be implemented to perform meshing of the device.
    %
    % Methods:
    %   Device(layer) - Constructor to create an instance of the Device class.
    %   setMaterial(varargin) - Sets the material for the device.
    %   combineDevice(varargin) - Combines multiple devices into one.
    %   changeUnit(unit) - Changes the unit of the device's properties.
    %
    % See also:
    %   Layout, Material, Grid, Axis, Grid2D, Device2D

    properties (SetAccess = protected)
        geo
        mat
        layer
        meshflag
        mesh
        eps
        mu
    end

    methods (Abstract)
        meshDevice
    end

    methods

        function obj = Device(layer)
            % Device constructor
            %   Constructs an instance of the Device class with a specified layer number.
            %
            %   Syntax:
            %     obj = Device(layer)
            %
            %   Input:
            %     layer - A positive integer specifying the layer number of the device.
            %
            %   Output:
            %     obj - An instance of the Device class (to be subclassed).

            arguments
                layer (1,1) {mustBeInteger, mustBePositive}
            end

            obj.geo = [];
            obj.mat = [];
            obj.layer(1) = layer;
            obj.meshflag = false;
            obj.mesh = [];
            obj.eps = [];
            obj.mu = [];
        end

        function setMaterial(obj, varargin)
            % setMaterial - Sets the material for the device
            %
            %   Syntax:
            %     setMaterial(varargin)
            %
            %   Input:
            %     varargin - Material properties or material name and wavelength.

            obj.mat = Material(varargin{:});
        end

        function obj = combineDevice(obj, varargin)
            % combineDevice - Combines multiple devices into one
            %
            %   Syntax:
            %     obj = combineDevice(varargin)
            %
            %   Input:
            %     varargin - Additional Device objects to be combined with the current object.
            %
            %   Output:
            %     obj - The combined Device object containing geometries, materials, and layers of all devices.

            arguments
                obj
            end
            arguments (Repeating)
                varargin
            end

            % check devices layers
            if numel(obj.layer) ~= numel(unique(obj.layer))
                dispError('Device:LayerOverlap');
            end

            % Combine devices into a single array
            device = [obj, varargin{:}];

            % Pre-allocate combined arrays
            combined_geo = {};
            combined_mat = [];
            combined_layer = [];

            % Iterate through each device
            for jj = 1:numel(device)
                % Add all geometries from the current device
                combined_geo = {combined_geo{:}, device(jj).geo};  % Concatenate cell arrays
                combined_mat = [combined_mat, device(jj).mat];
                combined_layer = [combined_layer, device(jj).layer];
            end

            obj.geo = combined_geo;
            obj.mat = combined_mat;
            obj.layer = combined_layer;
        end

    end

    methods (Static, Access = {?Device, ?FDFD})

        function s = YeeGridSampleVec(field_type, field_dim)
            % Yee 2x grid field components index
            vec_e = [1,1,1]; vec_h = [2,2,2];
            v1 = [1,0,0]; v2 = [0,1,0]; v3 = [0,0,1];
            vec_e1 = vec_e+v1; vec_e2 = vec_e+v2; vec_e3 = vec_e+v3;
            vec_h1 = vec_h-v1; vec_h2 = vec_h-v2; vec_h3 = vec_h-v3;

            % Select index vector
            var_name = sprintf('vec_%s%d', field_type, field_dim);
            s = eval(var_name);
        end

        function Ai = sampleMap(A_ds, field_type, dimension_idx)
            s = Device.YeeGridSampleVec(field_type, dimension_idx);

            switch numel(size(A_ds))
                case 1
                    Ai = A_ds(s(2):2:end);
                case 2
                    Ai = A_ds(s(1):2:end, s(2):2:end);
                case 3
                    Ai = A_ds(s(1):2:end, s(2):2:end, s(3):2:end);
            end
        end

        function [eps1, eps2, eps3, mu1, mu2, mu3] = sampleMapAll(eps_ds, mu_ds)
            eps1 = Device2D.sampleMap(eps_ds, 'e', 1); 
            eps2 = Device2D.sampleMap(eps_ds, 'e', 2); 
            eps3 = Device2D.sampleMap(eps_ds, 'e', 3);  

            mu1 = Device2D.sampleMap(mu_ds, 'h', 1); 
            mu2 = Device2D.sampleMap(mu_ds, 'h', 2); 
            mu3 = Device2D.sampleMap(mu_ds, 'h', 3);
        end

    end

    methods (Access = ?FDFD)

        function obj = changeUnit(obj, unit)
            arguments
                obj 
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})}
            end

            obj.mesh = obj.mesh.changeUnit(unit);
        end

    end

end