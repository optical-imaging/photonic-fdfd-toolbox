classdef Source1D < Source
    % Source1D: 1D FDFD input source for 2D simulation
    %
    % Description:
    %   The `Source1D` class defines a 1D input source for Finite Difference Frequency Domain (FDFD)
    %   simulations in 2D. It extends the base `Source` class to support various types of input sources
    %   such as eigenmodes, plane waves, and Gaussian beams.
    %
    % Properties:
    %   setflag - Boolean indicating if the source setup is complete.
    %   port - The port where the source is injected (Port1D object).
    %   mesh - The defined mesh grid for the source.
    %   lam - Wavelength of the source.
    %   neff - Effective refractive index of the source.
    %   Ex, Ey, Ez - Electric field components along the x, y, and z dimensions.
    %   Hx, Hy, Hz - Magnetic field components along the x, y, and z dimensions.
    %   phi - Injection phase of the source.
    %   theta - Rotation angle of the source in degrees.
    %   pol - Polarization of the source (e.g., TE, TM).
    %   label - Plane label indicating the simulation plane (e.g., 'xy', 'yz', 'zx').
    %
    % Methods:
    %   Source1D(plane_label, port, FDFD_mesh) - Constructor to create an instance of the Source1D class.
    %   setPhase(phi) - Sets the phase of the source.
    %   setAngle(theta_in_deg) - Sets the rotation angle of the source in degrees.
    %   printInfo() - Displays information about the source.
    %   assignEigenMode(eigenmode) - Assigns an existing eigenmode as the source.
    %   setEigenMode(device, wavelength, ...) - Solves an eigenmode and sets it as the source.
    %   setPlaneWave() - Sets a plane wave as the source (to be implemented).
    %   setGaussianBeam(wavelength, material_clad, w0, amplitude) - Sets a Gaussian beam as the source.
    %
    % See also:
    %   Source, Port1D, Device2D, Axis, Grid2D, Scattering

    properties (SetAccess = private)
        label
        pol
    end

    methods

        % Constructor
        function obj = Source1D(plane_label, port, FDFD_mesh)
            % Source1D constructor
            %   Constructs an instance of the Source1D class for a 2D simulation.
            %
            %   Syntax:
            %     obj = Source1D(plane_label, port, FDFD_mesh)
            %
            %   Input:
            %     plane_label - The label indicating the simulation plane (must be 'xy', 'yz', or 'zx').
            %     port - The port where the source is injected (Port1D object).
            %     FDFD_mesh - The mesh grid for the FDFD simulation (Grid2D object).
            %
            %   Output:
            %     obj - An instance of the Source1D class.

            arguments
                plane_label {mustBeMember(plane_label,{'xy', 'yz', 'zx'})}
                port Port1D
                FDFD_mesh Grid2D
            end

            % check port and plane
            if ~contains(plane_label, port.dir(2))
                dispError('Source1D:PortLableNotMatch', plane_label(1), plane_label(2));
            end

            obj@Source(port);
            obj.label = plane_label;
            obj.mesh = obj.getPortMesh(port, FDFD_mesh);
            obj.Ex = zeros(obj.mesh.n, 1);
            obj.Ey = zeros(obj.mesh.n, 1);
            obj.Ez = zeros(obj.mesh.n, 1);
            obj.Hx = zeros(obj.mesh.n, 1);
            obj.Hy = zeros(obj.mesh.n, 1);
            obj.Hz = zeros(obj.mesh.n, 1);
            obj.pol = [];
        end

        % Display
        function printInfo(obj)
            if obj.setflag
                info_name = {'Object Name';
                    'Wavelength';
                    'Phase';
                    'Polarization mode';
                    'Incidence direction'};
                value = {inputname(1);
                    dispWavelength(obj.lam);
                    num2str(obj.phi);
                    obj.pol;
                    obj.dir};

                maxNameLength = max(cellfun(@length, info_name));
                for ii = 1:numel(info_name)
                    fprintf('%*s: %-10s\n', maxNameLength+1, info_name{ii}, value{ii});
                end
            else
                dispMessage('Source1D:Initialization');
            end
        end

        function dispProfile(obj, field_name)
            arguments
                obj 
                field_name {mustBeMember(field_name, {'Ex','Ey','Ez','Hx','Hy','Hz'})}
            end

            field = obj.(field_name);

            figure;
            hold on; 
            legends = cell(1, size(field, 2));
            for ii = 1:size(field,2)
                plot(obj.mesh.v, field(:,ii), 'LineWidth', 1.5);
                legends{ii} = ['Source ' num2str(ii)];
            end
            hold off;
            xlabel('Port'), ylabel(field_name);
            legend(legends, 'Location', 'best');
        end

        % Other methods
        function obj = setEigenMode(obj, device, wavelength, varargin)
             % setEigenMode - Solves an eigenmode and sets it as the source
            %
            %   Syntax:
            %     obj = setEigenMode(device, wavelength, ...)
            %
            %   Input:
            %     device - The device for which to solve the eigenmode (Device2D object).
            %     wavelength - Wavelength of the source (positive real number).
            %     varargin - Additional arguments for solving the eigenmode.

            arguments
                obj
                device Device2D
                wavelength (1,1) double {mustBePositive}
            end

            arguments (Repeating)
                varargin
            end

            eigenmode = device.solveEigenModeSource(obj.port, wavelength, varargin{:});
            obj = assignEigenMode(obj,eigenmode);
            dispMessage('Source1D:VariableSetup');
        end

        function obj = setPlaneWave(obj, wavelength, material_clad, pol, amplitude)
            % setPlaneWave - Placeholder method to set a plane wave as the source
            
            arguments
                obj 
                wavelength (1,1) double {mustBePositive}
                material_clad (1,1) Material
                pol {mustBeMember(pol, {'TE', 'TM'})}
                amplitude (1, 1) {mustBePositive} = 1
            end

            % profile (constant)
            profile = amplitude*ones(size(obj.mesh.v'));

            % assign field
            field_pol = setdiff('xyz', obj.label);
            switch pol
                case 'TE'
                    field_name = ['E', field_pol];
                case 'TM'
                    field_name = ['H', field_pol];
            end
            obj.(field_name) = profile;

            obj.lam = wavelength;
            obj.pol = pol;
            obj.neff = obj.correctNeff(sqrt(material_clad.eps));
            obj.setflag = true;
            dispMessage('Source1D:VariableSetup');
        end

        function obj = setGaussianBeam(obj, wavelength, material_clad, w0, pol, amplitude)
            % setGaussianBeam - Sets a Gaussian beam as the source
            %
            %   Syntax:
            %     obj = setGaussianBeam(wavelength, material_clad, w0, amplitude)
            %
            %   Input:
            %     wavelength - Wavelength of the source (positive real number).
            %     material_clad - Cladding material (Material object).
            %     w0 - Beam waist of the Gaussian beam (positive real number).
            %     amplitude - Amplitude of the beam (positive real number, default is 1).
            
            arguments
                obj 
                wavelength (1,1) double {mustBePositive}
                material_clad (1,1) Material
                w0 (1,1) {mustBePositive}
                pol {mustBeMember(pol, {'TE', 'TM'})}
                amplitude (1, 1) {mustBePositive} = 1
            end

            % Gaussian beam profile
            port_center = (obj.mesh.v(1)+obj.mesh.v(end))/2;
            profile = amplitude*exp(-(obj.mesh.v'-port_center).^2/w0^2);

            % assign field
            field_pol = setdiff('xyz', obj.label);
            switch pol
                case 'TE'
                    field_name = ['E', field_pol];
                case 'TM'
                    field_name = ['H', field_pol];
            end
            obj.(field_name) = profile;

            obj.lam = wavelength;
            obj.pol = pol;
            obj.neff = obj.correctNeff(sqrt(material_clad.eps));
            obj.setflag = true;
            dispMessage('Source1D:VariableSetup');
        end


    end

    methods (Access = private)

        function obj = assignEigenMode(obj, eigenmode)
            arguments
                obj
                eigenmode EigenMode1D
            end

            obj.setflag = true;
            obj.lam = eigenmode.lam;
            obj.pol = eigenmode.pol;
            obj.neff = obj.correctNeff(eigenmode.neff);
            obj.Ex = eigenmode.Ex;
            obj.Ey = eigenmode.Ey;
            obj.Ez = eigenmode.Ez;
            obj.Hx = eigenmode.Hx;
            obj.Hy = eigenmode.Hy;
            obj.Hz = eigenmode.Hz;
        end

        function neff_sign = correctNeff(obj,neff)
            switch obj.dir(1)
                case '+'
                    neff_sign = neff;
                case '-'
                    neff_sign = -neff;
            end
        end

        function port_mesh = getPortMesh(obj, port, mesh)
            ind_axis_port = find(obj.label==obj.port.dir(2));
            switch ind_axis_port
                case 1
                    [~, ind_start] = min(abs(mesh.axis2.v-port.p1(2)));
                    [~, ind_end] = min(abs(mesh.axis2.v-port.p2(2)));
                    port_mesh = mesh.axis2.cutAxis(ind_start:ind_end);
                case 2
                    [~, ind_start] = min(abs(mesh.axis1.v-port.p1(1)));
                    [~, ind_end] = min(abs(mesh.axis1.v-port.p2(1)));
                    port_mesh = mesh.axis1.cutAxis(ind_start:ind_end);
            end
        end

    end

    methods (Access = ?FDFD2D)

        function obj = extendField(obj, axis)
            [~, ind_min] = min(abs(obj.mesh.v(1)-axis.v));
            [~, ind_max] = min(abs(obj.mesh.v(end)-axis.v));
            n_source = size(obj.Ex,2);

            field_list = {'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'};
            for ii = 1:numel(field_list)
                F = zeros(axis.n,n_source);
                F(ind_min:ind_max,:) = obj.(field_list{ii});
                obj.(field_list{ii}) = F;
            end
        end

    end

end