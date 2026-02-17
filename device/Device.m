classdef (Abstract, Hidden) Device < handle
    %Device: Abstract base class for simulation devices (geometry + material + meshed maps)
    %
    % Key Properties (SetAccess = protected):
    %   seq        - Device sequence / ID
    %   geo        - Geometry specification (Layout objects or containers)
    %   mat        - Material object(s) associated with geometry
    %   mesh       - Mesh object (Axis / Grid2D / Grid3D)
    %   eps        - Relative permittivity map sampled on Yee grid (stored with 3 comps)
    %   mu         - Relative permeability map sampled on Yee grid (stored with 3 comps)
    %
    % Key Properties (Dependent):
    %   meshflag   - True if mesh has been assigned (i.e., ~isempty(mesh))
    %
    % Key Methods (Abstract):
    %   setGeometry(...)        - Assign geometry into geo (subclass-defined)
    %   meshDevice(...)         - Rasterize geometry/material onto provided mesh (subclass-defined)
    %
    % Key Methods (Public):
    %   Device(seq)             - Construct empty device container with given seq
    %   setMaterial(...)        - Create and assign Material(...) to mat
    %   shiftMesh(dim,mode,val) - Shift internal mesh using Grid/Axes shiftGrid(...)
    %   setBitmap(eps,mu,ds,...) - Directly assign eps/mu bitmaps and construct mesh
    %
    % Key Methods (Static, Access = {?Device, ?FDFD)}):
    %   YeeGridSampleVec(t,dim) - Return Yee-grid sampling offset vector for E/H component
    %   sampleMap(A2x,t,dim)    - Sample a 2x-doublesampled map down to Yee locations
    %   sampleMapAll(eps2x,mu2x)- Sample eps/mu into (eps1,eps2,eps3,mu1,mu2,mu3)
    %
    % Key Methods (Access = ?Model):
    %   updateMesh(mesh)        - Replace internal mesh (workflow-controlled by Model)

    properties (SetAccess = protected)
        seq
        geo
        mat
        mesh
        eps
        mu
    end

    properties (Dependent)
        meshflag
    end

    methods (Abstract)
        setGeometry
        meshDevice
        setBitmap
    end

    methods
        % constructor
        function obj = Device(device_seq)
            arguments
                device_seq (1,1) {mustBeInteger, mustBePositive}
            end
            obj.geo = [];
            obj.mat = [];
            obj.seq = device_seq;
            obj.mesh = [];
            obj.eps = [];
            obj.mu = [];
        end

        % dependent
        function val = get.meshflag(obj)
            if ~isempty(obj.mesh)
                val = true;
            else
                val = false;
            end
        end

        % device set up
        function setMaterial(obj, varargin)
            obj.mat = Material(varargin{:});
        end

        % device meshing
        function obj = shiftMesh(obj, dim, mode, new_value)
            obj.mesh = obj.mesh.shiftGrid(dim, mode, new_value);
        end
    end

    methods (Static)
        function combined_device = combineDevice(varargin)
            arguments (Repeating)
                varargin
            end

            % Collect and validate inputs
            if isempty(varargin)
                dispError('Device:EmptyInput');
                return
            end

            device_list = [varargin{:}];

            % Type check (robust if users accidentally pass non-Device)
            if ~all(arrayfun(@(d) isa(d, 'Device'), device_list))
                dispError('Device:InvalidInputType');
            end

            % Sequence check & rearrange
            all_seq = [device_list.seq];
            if numel(all_seq) ~= numel(unique(all_seq))
                dispError('Device:SeqOverlap');
            end
            [~, order] = sort(all_seq, 'ascend');
            device_list = device_list(order);

            % Class consistency check (recommended)
            cls0 = class(device_list(1));
            if ~all(arrayfun(@(d) strcmp(class(d), cls0), device_list))
                dispError('Device:ClassMismatch');
            end

            % Create new combined device (same concrete class as inputs)
            combined_device = feval(cls0, 1);  % dummy seq; overwrite below

            % Aggregate
            combined_seq = [];
            combined_geo = {};
            combined_mat = [];

            for jj = 1:numel(device_list)
                combined_seq = [combined_seq, device_list(jj).seq];
                combined_geo = {combined_geo{:}, device_list(jj).geo};
                combined_mat = [combined_mat, device_list(jj).mat];
            end

            % Assign fields (allowed: within Device methods, protected set is permitted)
            combined_device.seq = combined_seq;
            combined_device.geo = combined_geo;
            combined_device.mat = combined_mat;

            % New composite should be re-meshed
            combined_device.mesh = [];
            combined_device.eps  = [];
            combined_device.mu   = [];
        end
    end

    methods (Static, Access = {?Device, ?FDFD})
        function s = YeeGridSampleVec(field_type, field_dim)
            vec_e = [1,1,1]; vec_h = [2,2,2];
            v1 = [1,0,0]; v2 = [0,1,0]; v3 = [0,0,1];
            vec_e1 = vec_e + v1; vec_e2 = vec_e + v2; vec_e3 = vec_e + v3;
            vec_h1 = vec_h - v1; vec_h2 = vec_h - v2; vec_h3 = vec_h - v3;
            var_name = sprintf('vec_%s%d', field_type, field_dim);
            s = eval(var_name);
        end

        function Ai = sampleMap(A_ds, field_type, dimension_idx)
            s = Device.YeeGridSampleVec(field_type, dimension_idx);
            switch ndims(A_ds)
                case 1
                    Ai = A_ds(s(1):2:end);
                case 2
                    Ai = A_ds(s(1):2:end, s(2):2:end);
                case 3
                    Ai = A_ds(s(1):2:end, s(2):2:end, s(3):2:end);
            end
        end

        function [eps1, eps2, eps3, mu1, mu2, mu3] = sampleMapAll(eps_ds, mu_ds)
            eps1 = Device.sampleMap(eps_ds, 'e', 1);
            eps2 = Device.sampleMap(eps_ds, 'e', 2);
            eps3 = Device.sampleMap(eps_ds, 'e', 3);
            mu1  = Device.sampleMap(mu_ds,  'h', 1);
            mu2  = Device.sampleMap(mu_ds,  'h', 2);
            mu3  = Device.sampleMap(mu_ds,  'h', 3);
        end
    end

    methods (Access = ?Model) % for solving workflow
        function updateMesh(obj,new_mesh)
            obj.mesh = new_mesh;
        end
    end
end
