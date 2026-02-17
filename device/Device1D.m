classdef Device1D < Device
    %Device1D: 1D simulation device (1D geometry + 1D meshing into eps/mu maps)
    %
    % Key Properties (inherited, 1D meaning):
    %   seq        - Device ID
    %   geo        - 1D geometry specification (subclass-defined)
    %   mat        - Material definition for 1D regions
    %   mesh       - Axis object for the 1D coordinate
    %   eps        - Permittivity map (expanded internally to include 3 components)
    %   mu         - Permeability map (expanded internally to include 3 components)
    %   meshflag   - True when mesh is assigned
    %
    % Key Methods (implements abstract API):
    %   setGeometry(...)        - Assign 1D geometry definition into geo
    %   meshDevice(axis)        - Rasterize 1D structure onto an Axis mesh
    %
    % Key Methods (inherited utility):
    %   setMaterial(...)        - Assign Material to this device
    %   shiftMesh(... )         - Shift Axis by min/max/center alignment
    %   setBitmap(...)          - Directly define eps/mu and build an Axis mesh

    methods
        % constructor
        function obj = Device1D(device_seq)
            arguments
                device_seq (1,1) {mustBeInteger, mustBePositive}
            end
            obj@Device(device_seq);
        end

        % device set up
        function obj = setGeometry(obj, start_end_point)
            arguments
                obj
                start_end_point (1,2) double
            end
            p1 = start_end_point(1);
            p2 = start_end_point(2);
            if p1 >= p2
                dispError('Device1D:InvalidGeometry');
            end
            obj.geo = Line1D(p1, p2);
        end

        % device meshing
        function meshDevice(obj, axis, bg_mat)
            arguments
                obj
                axis Axis
                bg_mat Material = []
            end

            if isempty(bg_mat)
                eps_bg = 1;
                mu_bg = 1;
            else
                eps_bg = bg_mat.eps;
                mu_bg = bg_mat.mu;
            end

            if isempty(obj.geo)
                geolist = {};
            elseif iscell(obj.geo)
                geolist = obj.geo;
            elseif isvector(obj.geo)
                geolist = num2cell(obj.geo);
            else
                geolist = {obj.geo};
            end

            % 2x grid
            axis_ds = axis.doublesampleAxis;
            eps_ds = eps_bg*ones(axis_ds.n, 1);
            mu_ds  = mu_bg*ones(axis_ds.n, 1);

            for ii = 1:numel(geolist)
                g = geolist{ii};
                in = axis_ds.v >= g.p1(1) & axis_ds.v <= g.p2(1);
                eps_val = 1;
                mu_val  = 1;
                if ii <= numel(obj.mat)
                    m = obj.mat(ii);
                    if isnumeric(m)
                        eps_val = real(m);
                        mu_val  = 1;
                    else
                        if isprop(m,'eps') && ~isempty(m.eps), eps_val = m.eps; end
                        if isprop(m,'mu')  && ~isempty(m.mu),  mu_val  = m.mu;  end
                    end
                end
                eps_ds(in) = eps_val;
                mu_ds(in) = mu_val;
            end

            [eps1, eps2, eps3, mu1, mu2, mu3] = Device.sampleMapAll(eps_ds, mu_ds);

            obj.meshflag = true;
            obj.mesh = axis;
            obj.eps = cat(2,eps1,eps2,eps3);
            obj.mu  = cat(2,mu1,mu2,mu3);
        end

        function setBitmap(obj, eps_bit, mu_bit, step_size, axis_label)
            arguments
                obj
                eps_bit (:,1) double {mustBeRealAndGE1}
                mu_bit (:,1) double {mustBeRealAndGE1}
                step_size (1,1) double {mustBeReal, mustBePositive}
                axis_label char {mustBeMember(axis_label, {'x','y','z'})}
            end
            % eps and mu size match
            if ~isequal(size(eps_bit), size(mu_bit))
                dispError('Device:BitmapSizeMismatch');
            end

            obj.eps = repmat(eps_bit, [1 3]);
            obj.mu  = repmat(mu_bit, [1 3]);
            obj.mesh = Axis(size(eps_bit,1), step_size, axis_label);
        end
    end

    methods (Access = ?Source1D) % for port-cut device construction
        function setTotalBitmap(obj, eps_bit, mu_bit, axis)
            arguments
                obj
                eps_bit (:,3) double {mustBeRealAndGE1}
                mu_bit (:,3) double {mustBeRealAndGE1}
                axis Axis
            end
            obj.eps = eps_bit;
            obj.mu = mu_bit;
            obj.mesh = axis;
        end
    end
end
