classdef Device2D < Device
    % Device2D: Derived class for 2D devices in the simulation.
    %
    % Description:
    %   The `Device2D` class is a subclass of the abstract `Device` class, providing implementation
    %   for devices specifically defined in a two-dimensional space. It includes functionalities
    %   to set geometry, solve eigenmodes, and compute scattering parameters.
    %
    % Properties:
    %   label - Label indicating the plane of the 2D device ('xy', 'yz', or 'zx').
    %
    % Methods:
    %   Device2D(layer, plane_label) - Constructor to create an instance of the Device2D class.
    %   setGeometry(shape_name, varargin) - Sets the geometry of the device.
    %   solveEigenMode(wavelength, num_mode) - Solves the eigenmode for a given wavelength and number of modes.
    %   solveEigenModeSource(port, wavelength, varargin) - Solves the eigenmode for a given source port and wavelength.
    %   solveScattering(source, pml_thickness) - Solves the scattering parameters for the device.
    %   combineDevice(varargin) - Combines multiple 2D devices into one.
    %   meshDevice(grid) - Performs meshing on the 2D device using the given grid.
    %   dispImg(plot_part) - Displays the geometry of the device with the specified material property.
    %
    % See also:
    %   Device, Layout2D, Material, Grid, Axis, Grid2D, Port1D

    properties
        label
    end

    methods

        function obj = Device2D(layer, plane_label)
            % Device2D constructor
            %   Constructs an instance of the Device2D class with a specified layer and plane label.
            %
            %   Syntax:
            %     obj = Device2D(layer, plane_label)
            %
            %   Input:
            %     layer - A positive integer specifying the layer number of the device.
            %     plane_label - A string specifying the plane of the device ('xy', 'yz', 'zx').
            %
            %   Output:
            %     obj - An instance of the Device2D class.

            arguments
                layer (1,1) {mustBeInteger, mustBePositive}
                plane_label {mustBeMember(plane_label, {'xy', 'yz', 'zx'})}
            end

            obj@Device(layer);
            obj.label = plane_label;
        end

        % Set properties
        function setGeometry(obj, shape_name, varargin)
            % setGeometry - Sets the geometry of the device
            %
            %   Syntax:
            %     setGeometry(shape_name, varargin)
            %
            %   Input:
            %     shape_name - A string specifying the shape of the geometry ('Rectangle', 'Disk', 'Polygon', 'Ring').
            %     varargin - Additional parameters for the geometry (e.g., radius of disk).

            arguments
                obj
                shape_name {mustBeMember(shape_name, {'Rectangle', 'Disk', 'Polygon', 'Ring', 'GDSII'})}
            end
            arguments (Repeating)
                varargin
            end

            switch shape_name
                case 'Rectangle'
                    obj.geo = Rectangle(varargin{:});
                case 'Disk'
                    obj.geo = Disk(varargin{:});
                case 'Polygon'
                    obj.geo = Polygon(varargin{:});
                case 'Ring'
                    obj.geo = Ring(varargin{:});
                case 'GDSII'
                    % stack = dbstack('-completenames');
                    % callerpath = fileparts(stack(2).file);
                    % gdsrelativepath = varargin{:};  
                    % gdsfilepath = fullfile(callerpath, gdsrelativepath);
                    obj.geo = GDSII(varargin{:});
            end
        end

        % FDFD simulation
        function eigenmode = solveEigenMode(obj, wavelength, num_mode)
            % solveEigenMode - Solves the eigenmode for the 2D device
            %
            %   Syntax:
            %     eigenmode = solveEigenMode(wavelength, num_mode)
            %
            %   Input:
            %     wavelength - Wavelength for which to solve the eigenmode.
            %     num_mode - Number of mode to solve for.
            %
            %   Output:
            %     eigenmode - The solved eigenmode.

            if ~obj.meshflag
                dispError('Device2D:MeshDevice');
            end

            eigenmode = EigenMode2D(obj.eps, obj.mu, obj.mesh, obj.label, wavelength, num_mode);
        end

        function eigenmode = solveEigenModeSource(obj, port, wavelength, varargin)
            % solveEigenModeSource - Solves the eigenmode a given source port on the 2D device
            %
            %   Syntax:
            %     eigenmode = solveEigenModeSource(port, wavelength, varargin)
            %
            %   Input:
            %     port - The source port (Port1D) for which to solve the eigenmode.
            %     wavelength - Wavelength for which to solve the eigenmode.
            %     varargin - Additional parameters for the eigenmode solver.
            %
            %   Output:
            %     eigenmode - The solved 1D eigenmode for source.

            arguments
                obj
                port Port1D
                wavelength (1,1) double {mustBeReal, mustBePositive}
            end
            arguments (Repeating)
                varargin
            end

            if ~obj.meshflag
                dispError('Device2D:MeshDevice');
            end

            [eps_port, mu_port, axis_port, label_port] = obj.CrossSection(port);
            eigenmode = EigenMode1D(eps_port, mu_port, axis_port, label_port, wavelength, varargin{:});
        end

        function scattering = solveScattering(obj, source, pml_thickness)
            % solveScattering - Solves the scattering parameters for the device
            %
            %   Syntax:
            %     scattering = solveScattering(source, pml_thickness)
            %
            %   Input:
            %     source - The source object (Source1D) for the scattering simulation.
            %     pml_thickness - Thickness of the PML (Perfectly Matched Layer) in grid points.
            %
            %   Output:
            %     scattering - The solved scattering results.

            arguments
                obj
                source Source1D
                pml_thickness (1,2) {mustBePositive, mustBeInteger}
            end

            scattering = Scattering2D(obj, source, pml_thickness);
        end

        % Other methods
        function obj = combineDevice(obj, varargin)
            % combineDevice - Combines multiple 2D devices into one
            %
            %   Syntax:
            %     obj = combineDevice(varargin)
            %
            %   Input:
            %     varargin - Additional Device2D objects to be combined with the current object.
            %
            %   Output:
            %     obj - The combined Device2D object containing geometries, materials, and layers of all devices.

            arguments
                obj
            end
            arguments (Repeating)
                varargin Device2D
            end

            obj = combineDevice@Device(obj, varargin{:});
        end

        function obj = meshDevice(obj, grid)
            % meshDevice - Performs meshing on the 2D device using the given grid
            %
            %   Syntax:
            %     obj = meshDevice(grid)
            %
            %   Input:
            %     grid - The grid object (Grid2D) used to perform the meshing.

            arguments
                obj
                grid Grid2D
            end

            % reconstruct 2x axis
            axis1_ds = grid.axis1.doublesampleAxis;
            axis2_ds = grid.axis2.doublesampleAxis;
            [AXIS1, AXIS2] = ndgrid(axis1_ds.v, axis2_ds.v);

            % initializa eps and mu map
            eps_ds = ones(axis1_ds.n, axis2_ds.n);
            mu_ds = ones(axis1_ds.n, axis2_ds.n);

            for ii = 1:numel(obj.geo)
                geometry_map = inpolygon(AXIS1, AXIS2,...
                    obj.geo{ii}.s.Vertices(:, 1), obj.geo{ii}.s.Vertices(:, 2));

                eps_ds(geometry_map) = obj.mat(ii).eps;
                mu_ds(geometry_map) = obj.mat(ii).mu;
            end
            [eps1, eps2, eps3, mu1, mu2, mu3] = Device.sampleMapAll(eps_ds, mu_ds);
            eps = cat(3, eps1, eps2, eps3);
            mu = cat(3, mu1, mu2, mu3);

            obj.meshflag = true;
            obj.mesh = grid;
            obj.eps = eps;
            obj.mu = mu;
        end

        % Display
        function dispImg(obj, plot_part)
            arguments
                obj
                plot_part {mustBeMember(plot_part, {'eps', 'mu'})} = 'eps'
            end
            % Define color range
            minValue = 1;
            maxValue = 15;
            cmap = jet(256); % Colormap from blue to red

            figure;
            hold on;
            axis image;

            switch plot_part
                case 'eps'
                    title([inputname(1), ' \epsilon_r']);
                case 'mu'
                    title([inputname(1), ' \mu_r']);
            end

            for ii = 1:numel(obj.geo)
                % Normalize material value to get color index
                switch plot_part
                    case 'eps'
                        matValue = real(obj.mat(ii).eps);
                    case 'mu'
                        matValue = real(obj.mat(ii).mu);
                end
                colorIdx = round((matValue - minValue) / (maxValue - minValue) * 255) + 1;
                colorIdx = max(1, min(colorIdx, 256)); % Clamp to valid range

                plot(obj.geo{ii}.s, 'FaceColor', cmap(colorIdx, :), 'FaceAlpha', 1, 'EdgeColor', 'none');
                xlabel(obj.label(1)), ylabel(obj.label(2)), axis image;
            end

            hold off;
            colorbar;
            colormap(cmap);
            clim([minValue, maxValue]); % Set color limits to match material range
        end
    end

    methods (Access = private)

        function [eps_cs, mu_cs, axis_cs, label_cs]  = CrossSection(obj, port)
            arguments
                obj
                port Port1D
            end

            ind_axis_port = find(obj.label==port.dir(2));
            if isempty(ind_axis_port)
                dispError('Device2D:WrongPortDirection', obj.label(1), obj.label(2));
            end
            switch ind_axis_port
                case 1
                    [~, ind_start] = min(abs(obj.mesh.axis2.v-port.p1(2)));
                    [~, ind_end] = min(abs(obj.mesh.axis2.v-port.p2(2)));
                    [~, ind_position] = min(abs(obj.mesh.axis1.v-port.p1(1)));

                    eps_cs = obj.eps(ind_position, ind_start:ind_end, :);
                    mu_cs = obj.mu(ind_position, ind_start:ind_end, :);
                    axis_cs = obj.mesh.axis2.cutAxis(ind_start:ind_end);
                case 2
                    [~, ind_start] = min(abs(obj.mesh.axis1.v-port.p1(1)));
                    [~, ind_end] = min(abs(obj.mesh.axis1.v-port.p2(1)));
                    [~, ind_position] = min(abs(obj.mesh.axis2.v-port.p1(2)));

                    eps_cs = obj.eps(ind_start:ind_end, ind_position, :);
                    mu_cs = obj.mu(ind_start:ind_end, ind_position, :);
                    axis_cs = obj.mesh.axis1.cutAxis(ind_start:ind_end);
            end
            label_cs = setdiff(obj.label,port.dir(2));
        end

    end

    methods (Access = ?Model2D)

        function device_copy = copyDevice(obj)
            device_copy = Device2D(1, obj.label);
            device_copy.geo = obj.geo;
            device_copy.mat = obj.mat;
            device_copy.layer = obj.layer;
            device_copy.meshflag = obj.meshflag;
            device_copy.mesh = obj.mesh;
            device_copy.eps = obj.eps;
            device_copy.mu = obj.mu;
        end

    end


end