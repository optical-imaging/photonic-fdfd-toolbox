classdef Source1D < Source
    %Source1D: 1D port-defined excitation source for 2D FDFD scattering simulations
    %
    % Key Properties (inherited from Source):
    %   lam            - Source wavelength (m)
    %   port           - Port1D defining injection line and direction
    %   mesh           - 1D port mesh (Axis) extracted from a Grid2D
    %   neff           - Effective index of injected mode (sign corrected by dir)
    %   phi            - Phase offset (rad)
    %   theta          - Incident angle (deg) (stored; not used by current implementations here)
    %   Ex,Ey,Ez       - Electric fields on port mesh, size Np×Ns
    %   Hx,Hy,Hz       - Magnetic fields on port mesh, size Np×Ns
    %
    % Key Properties (Source1D-specific):
    %   pol            - Polarization: 'TE'|'TM'
    %
    % Dependent Properties:
    %   label          - Plane label forwarded from port.label (e.g., 'xy','yz','zx')
    %   dir            - Injection direction forwarded from port.dir (inherited dependent)
    %   setflag        - True if mesh/neff and all field components are populated (inherited dependent)
    %
    % Key Methods (inherited from Source, still used here):
    %   Source(lam,port)                 - Construct base Source container and defaults
    %   get.dir()                        - Return port.dir
    %   get.setflag()                    - Return readiness flag (fields + neff + mesh)
    %   setPhase(phi)                    - Set phase offset
    %   setAngle(theta_deg)              - Set incident angle (deg)
    %   setAmplitude(amp,seq,'E'|'H')    - Scale selected source instance(s) by field magnitude
    %
    % Key Methods (implements Source abstract API):
    %   Source1D(lam,port,pol)           - Construct 1D source with polarization
    %   get.label()                      - Return port.label
    %   setMesh(grid2d)                  - Extract 1D port Axis mesh from Grid2D (label must match)
    %   setEigenMode(device2d,...)       - Solve EigenMode1D on port cross-section and assign fields
    %   setPlaneWave(BG_material)        - Create uniform-profile plane-wave fields on port mesh
    %   setGaussianBeam(BG_material,w0)  - Create Gaussian-profile fields on port mesh
    %
    % Key Methods (display/debug):
    %   printInfo()                      - Print basic source settings if initialized
    %   dispProfile(field_name)          - Plot selected field component(s) vs port coordinate
    %
    % Key Methods (Static, protected; required by Source):
    %   getPortMesh(port,grid2d)         - Return Axis along the transverse direction of the port
    %
    % Key Methods (Access = {?Scattering, ?Source}; required by Source):
    %   extendField(axis)                - Return struct of fields extended onto a target Axis

    properties (SetAccess = private)
        pol
    end

    properties (Dependent)
        label
    end

    methods
        % conostructor
        function obj = Source1D(wavelength, port, pol)
            arguments
                wavelength (1,1) double {mustBePositive}
                port Port1D
                pol string {mustBeMember(pol, {'TE', 'TM'})}
            end

            obj@Source(wavelength, port);
            obj.pol = pol;
        end

        % dependent
        function val = get.label(obj)
            val = obj.port.label;
        end 

        % set up mesh based on port
        function setMesh(obj, mesh2d)
            arguments
                obj 
                mesh2d Grid2D
            end
            if ~strcmp(mesh2d.label, obj.label)
                dispError('Source1D:SourceLabelNotInMesh');
            end
            obj.mesh = Source1D.getPortMesh(obj.port, mesh2d);
        end

        % source field types
        function setEigenMode(obj, device, num_mode, eigenmode1d_flip)
            arguments
                obj
                device Device2D
                num_mode (1,:) double {mustBeInteger, mustBePositive} = 1
                eigenmode1d_flip (1,1) logical = false
            end

            % create port-cut device1D
            [eps_port, mu_port, ~] = device.CrossSection(obj.port);
            device_port = Device1D(1);
            device_port.setTotalBitmap(eps_port, mu_port, obj.mesh);

            % solve port-cut device eigenmode
            eigenmode_src = EigenMode1D(obj.lam, num_mode, obj.pol);
            eigenmode_src.setMesh(obj.mesh);
            eigenmode_src.solveFDFD(device_port, eigenmode1d_flip);

            obj.assignEigenMode(eigenmode_src);
        end

        function setPlaneWave(obj, BG_material)
            arguments
                obj
                BG_material (1,1) Material
            end

            % initialial field components
            obj.Ex = zeros(obj.mesh.n, 1); % for plane wave, there is only one mode
            obj.Ey = zeros(obj.mesh.n, 1);
            obj.Ez = zeros(obj.mesh.n, 1);
            obj.Hx = zeros(obj.mesh.n, 1);
            obj.Hy = zeros(obj.mesh.n, 1);
            obj.Hz = zeros(obj.mesh.n, 1);

            % assign fields
            profile = ones(obj.mesh.n,1);
            field_pol = setdiff('xyz', obj.label);
            switch obj.pol
                case 'TE'
                    field_name = ['E', field_pol];
                case 'TM'
                    field_name = ['H', field_pol];
            end
            obj.(field_name) = profile;

            obj.correctNeff(sqrt(BG_material.eps*BG_material.mu));
        end


        function setGaussianBeam(obj, BG_material, waist)
            arguments
                obj
                BG_material (1,1) Material
                waist (1,1) double {mustBeGreaterThan(waist, 0)}
            end

            % Gaussian pulse
            center = mean([obj.mesh.v(1), obj.mesh.v(end)]);
            profile = exp(-(obj.mesh.v' - center).^2/waist^2);

            % initialial field components
            obj.Ex = zeros(obj.mesh.n, 1);
            obj.Ey = zeros(obj.mesh.n, 1);
            obj.Ez = zeros(obj.mesh.n, 1);
            obj.Hx = zeros(obj.mesh.n, 1);
            obj.Hy = zeros(obj.mesh.n, 1);
            obj.Hz = zeros(obj.mesh.n, 1);

            % assign fields
            field_pol = setdiff('xyz', obj.label);
            switch obj.pol
                case 'TE'
                    field_name = ['E', field_pol];
                case 'TM'
                    field_name = ['H', field_pol];
            end
            obj.(field_name) = profile;

            obj.correctNeff(sqrt(BG_material.eps*BG_material.mu));
        end

        % display
        function printInfo(obj)
            if obj.setflag
                info_name = {'Object Name', 'Wavelength', 'Phase', ...
                    'Polarization mode', 'Incidence direction'};
                value = {inputname(1), dispWavelength(obj.lam), ...
                    num2str(obj.phi), obj.pol, obj.dir};

                maxlen = max(cellfun(@length, info_name));
                for ii = 1:numel(info_name)
                    fprintf('%*s: %-10s\n', maxlen + 1, info_name{ii}, value{ii});
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
            figure; hold on;
            for ii = 1:size(field, 2)
                plot(obj.mesh.v, field(:, ii), 'LineWidth', 1.5);
            end
            hold off;
            xlabel('Port'); ylabel(field_name);
            legend("Source " + (1:size(field,2)), 'Location', 'best');
        end

    end

    methods (Static, Access = protected)
        function port_mesh = getPortMesh(port, mesh)
            ind_axis_port = find(mesh.label == port.dir(2));
            switch ind_axis_port
                case 1
                    [~, i1] = min(abs(mesh.axis2.v - port.p1(2)));
                    [~, i2] = min(abs(mesh.axis2.v - port.p2(2)));
                    port_mesh = mesh.axis2.cutAxis(i1:i2);
                case 2
                    [~, i1] = min(abs(mesh.axis1.v - port.p1(1)));
                    [~, i2] = min(abs(mesh.axis1.v - port.p2(1)));
                    port_mesh = mesh.axis1.cutAxis(i1:i2);
            end
        end

    end

    methods (Access = {?Scattering, ?Source})
        function extended_field = extendField(obj, axis)
            [~, i1] = min(abs(obj.mesh.v(1) - axis.v));
            [~, i2] = min(abs(obj.mesh.v(end) - axis.v));
            nsrc = size(obj.Ex,2);

            extended_field = struct();
            extended_field.Ex = [];
            extended_field.Ey = [];
            extended_field.Ez = [];
            extended_field.Hx = [];
            extended_field.Hy = [];
            extended_field.Hz = [];
            for field = ["Ex","Ey","Ez","Hx","Hy","Hz"]
                F = zeros(axis.n, nsrc);
                F(i1:i2, :) = obj.(field);
                extended_field.(field) = F;
            end
        end
    end

end
