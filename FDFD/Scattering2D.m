classdef (Hidden) Scattering2D < FDFD2D
    % Scattering2D: 2D scattering FDFD simulation
    %
    % Description:
    %   The `Scattering2D` class is used for performing two-dimensional scattering simulations using finite-difference
    %   frequency-domain (FDFD) techniques. This class allows users to analyze how electromagnetic waves scatter
    %   from a given device in a 2D domain.
    %   This class inherits from `FDFD2D` and provides additional methods and properties specific to scattering analysis.
    %
    % Properties:
    %   mesh - The defined mesh grid for the 2D domain.
    %   lam - The wavelength of the electromagnetic wave in the simulation.
    %   solflag - A flag indicating if the solution has been computed.
    %   Ex, Ey, Ez - Electric field components along x, y, and z dimensions.
    %   Hx, Hy, Hz - Magnetic field components along x, y, and z dimensions.
    %   pol - Polarization of the wave ('TE' or 'TM').
    %
    % Methods:
    %   Scattering2D(device, source, PML_thickness) - Constructor to create a 2D scattering solver object.
    %
    % Example:
    %   % Define a device and source for scattering analysis
    %   device = Device2D(Grid2D(Axis('x', 0, 1, 100), Axis('y', 0, 1, 100)), eps, mu);
    %   source = Source1D(...); % Define a source appropriately
    %   PML_thickness = [10, 10];
    %
    %   % Create a Scattering2D solver object
    %   scat = Scattering2D(device, source, PML_thickness);
    %
    % Notes:
    %   - This class is marked as "Hidden", which means it is not intended to be directly accessed by users.
    %   - Ensure that the `Device2D` and `Source1D` objects are properly defined before using them.
    %
    % See Also:
    %   FDFD, FDFD2D, Device2D, Source1D

    properties
        pol
    end

    properties (Constant)
        scamp = 3;
        conduct = 300;
        power = 3;
    end

    methods

         function obj = Scattering2D(device, source, PML_thickness)
            % Scattering2D: Constructor to create a 2D scattering solver object
            %
            %   Syntax:
            %     obj = Scattering2D(device, source, PML_thickness)
            %
            %   Input:
            %     device - A `Device2D` object representing the scattering object.
            %     source - A `Source1D` object representing the electromagnetic source.
            %     PML_thickness - Thickness of the perfectly matched layers (PML) in grid points (2-element positive integer array).
            %
            %   Output:
            %     obj - An instance of the `Scattering2D` class.

            arguments
                device Device2D
                source Source1D
                PML_thickness (1,2) {mustBePositive, mustBeInteger}
            end

            % pre-check
            if ~device.meshflag
                dispError('Device2D:MeshDevice');
            end

            if ~source.setflag
                dispError('FDFD2D:DefineSourceFirst');
            end

            ind_axis_port = find(device.label==source.port.dir(2));
            switch ind_axis_port
                case 1
                    axis_source = device.mesh.axis2;
                    axis_prop = device.mesh.axis1;
                    mesh_tp_flag = false;
                case 2
                    axis_source = device.mesh.axis1;
                    axis_prop = device.mesh.axis2;
                    mesh_tp_flag = true;
            end

            obj@FDFD2D(device.mesh);

            % constant
            eta0 = Constant('eta0').v; % vacuum impedance
            k0 = 2*pi/source.lam; % wave number

            % change all units to m
            obj.mesh = device.mesh;
            device = device.changeUnit('m');

            % prepare source
            source = source.extendField(axis_source);
            field_pol = setdiff('xyz', device.label, 'stable');
            switch source.pol
                case 'TE'
                    source_field = source.(['E',field_pol]);
                case 'TM'
                    source_field = source.(['H',field_pol]);
            end

            % building FDFD matrices
            % 2x grid
            axis1_ds = device.mesh.axis1.doublesampleAxis;
            axis2_ds = device.mesh.axis2.doublesampleAxis;
            [AXIS1, AXIS2] = ndgrid(axis1_ds.v, axis2_ds.v);
            switch mesh_tp_flag
                case false
                    V_prop = AXIS1;
                    V_source = AXIS2;
                case true
                    V_prop = AXIS2;
                    V_source = AXIS1;
            end

            % derivative matrices
            [DE1, DE2] = FDFD.DM2D(device.mesh, 'E');
            [DH1, DH2] = FDFD.DM2D(device.mesh, 'M');
            DE1 = DE1/k0;
            DE2 = DE2/k0;
            DH1 = DH1/k0;
            DH2 = DH2/k0;

            % SCPML
            PML.scamp = obj.scamp;
            PML.conduct = obj.conduct;
            PML.power = obj.power;
            PML.thickness = PML_thickness;
            [S1E2_inv, S1E3_inv, S1H2_inv, S1H3_inv,...
             S2E1_inv, S2E3_inv, S2H1_inv, S2H3_inv] ...
                = Scattering2D.SCPML(2*device.mesh.axis1.n, 2*device.mesh.axis2.n, PML);

            % A matrix
            switch source.pol
                case 'TE'
                    mu1_inv = sparse(1./device.mu(:,:,1)); mu1_inv = diag(mu1_inv(:));
                    mu2_inv = sparse(1./device.mu(:,:,2)); mu2_inv = diag(mu2_inv(:));
                    eps3 = device.eps(:,:,3); eps3 = diag(sparse(eps3(:)));
                    A = S1E3_inv*DH1*mu2_inv*S1H2_inv*DE1 + S2E3_inv*DH2*mu1_inv*S2H1_inv*DE2 + eps3;
                case 'TM'
                    eps1_inv = sparse(1./device.eps(:,:,1)); eps1_inv = diag(eps1_inv(:));
                    eps2_inv = sparse(1./device.eps(:,:,2)); eps2_inv = diag(eps2_inv(:));
                    mu3 = device.mu(:,:,3); mu3 = diag(sparse(mu3(:)));
                    A = S1H3_inv*DE1*eps2_inv*S1E2_inv*DH1 + S2H3_inv*DE2*eps1_inv*S2E1_inv*DH2 + mu3;
            end


            % Q masking matrix
            Q = true(device.mesh.axis1.n, device.mesh.axis2.n);
            switch mesh_tp_flag
                case false
                    [~, ind_p] = min(abs(axis_prop.v-source.port.p1(1)));
                    Q(ind_p:end,:) = false;
                case ture
                    [~, ind_p] = min(abs(axis_prop.v-source.port.p1(2)));
                    Q(:,ind_p:end) = false;
            end

            switch source.dir(1)
                case '+'
                    Q = diag(sparse(Q(:)));
                case '-'
                    Q = diag(sparse(~Q(:)));
            end

            % incorporate source field
            n_source = numel(source.neff);
            b = zeros(device.mesh.axis1.n*device.mesh.axis2.n,n_source);
            for ii = 1:n_source
                fsrc_profile = FDFD2D.rotateProfile(source_field(:,ii), axis_source, axis_prop, source.theta);
                if mesh_tp_flag
                    fsrc_profile = fsrc_profile.';
                end

                fsrc_phase_ds = FDFD2D.rotatePhase(source.neff(ii)*k0,...
                    V_source, V_prop, axis_prop.v(ind_p), source.theta, source.phi);
                switch source.pol
                    case 'TE'
                        fsrc = fsrc_profile.*fsrc_phase_ds(1:2:end, 1:2:end);
                    case 'TM'
                        fsrc = fsrc_profile.*fsrc_phase_ds(2:2:end, 2:2:end);
                end

                b(:,ii) = (Q*A-A*Q)*fsrc(:);
            end

            % solve field
            f = zeros(device.mesh.axis1.n*device.mesh.axis2.n, n_source);
            for kk = 1:n_source
                b_sub = b(:,kk);
                f(:,kk) = A\b_sub;
            end

            % assign field value
            e1 = zeros(device.mesh.axis1.n*device.mesh.axis2.n, n_source);
            e2 = zeros(device.mesh.axis1.n*device.mesh.axis2.n, n_source);
            e3 = zeros(device.mesh.axis1.n*device.mesh.axis2.n, n_source);
            h1 = zeros(device.mesh.axis1.n*device.mesh.axis2.n, n_source);
            h2 = zeros(device.mesh.axis1.n*device.mesh.axis2.n, n_source);
            h3 = zeros(device.mesh.axis1.n*device.mesh.axis2.n, n_source);

            switch source.pol
                case 'TE'
                    e3 = f;
                    h1 = mu1_inv*S2H1_inv*DE2*e3;
                    h2 = -mu2_inv*S1H2_inv*DE1*e3;
                case 'TM'
                    h3 = f;
                    e1 = eps1_inv*S2E1_inv*DH2*h3;
                    e2 = -eps2_inv*S1E2_inv*DH1*h3;
            end

            % assign fields
            E_field_list = {'Ex', 'Ey', 'Ez'};
            H_field_list = {'Hx', 'Hy', 'Hz'};
            field_seq = Scattering2D.setFieldSeq(device.label);
            obj.(E_field_list{field_seq(1)}) = reshape(e1,...
                [device.mesh.axis1.n, device.mesh.axis2.n, n_source]);
            obj.(H_field_list{field_seq(1)}) = reshape(h1,...
                [device.mesh.axis1.n, device.mesh.axis2.n, n_source])/(-1i*eta0);
            obj.(E_field_list{field_seq(2)}) = reshape(e2,...
                [device.mesh.axis1.n, device.mesh.axis2.n, n_source]);
            obj.(H_field_list{field_seq(2)}) = reshape(h2,...
                [device.mesh.axis1.n, device.mesh.axis2.n, n_source])/(-1i*eta0);
            obj.(E_field_list{field_seq(3)}) = reshape(e3,...
                [device.mesh.axis1.n, device.mesh.axis2.n, n_source]);
            obj.(H_field_list{field_seq(3)}) = reshape(h3,...
                [device.mesh.axis1.n, device.mesh.axis2.n, n_source])/(-1i*eta0);

            obj.solflag = true;
            obj.lam = source.lam;
            obj.pol = source.pol;
         end

    end

    methods (Static, Access = private)

        function [S1E2_inv, S1E3_inv, S1H2_inv, S1H3_inv,...
                  S2E1_inv, S2E3_inv, S2H1_inv, S2H3_inv] = SCPML(N1_ds, N2_ds, PML)
            % PML settings
            scale_amp = PML.scamp;
            boundary_conduct = PML.conduct;
            power = PML.power;
            thickness = PML.thickness;

            % build SCPML layers
            thickness_ds = 2*thickness; % solve it in 2x grid

            s1_ds = ones(N1_ds, N2_ds);
            s2_ds = ones(N1_ds, N2_ds);

            n1 = (1:thickness_ds(1))';
            NN1 = repmat(n1,1,N2_ds);
            a1 = 1+scale_amp*(NN1/thickness_ds(1)).^power;
            c1 = boundary_conduct*sin(0.5*pi*NN1/thickness_ds(1)).^2;
            s1_ds(thickness_ds(1)-n1+1,:) = a1.*(1-1i*c1);
            s1_ds(N1_ds-thickness_ds(1)+n1, :) = a1.*(1-1i*c1);

            n2 = 1:thickness_ds(2);
            NN2 = repmat(n2,N1_ds,1);
            a2 = 1+scale_amp*(NN2/thickness_ds(2)).^power;
            c2 = boundary_conduct*sin(0.5*pi*NN2/thickness_ds(2)).^2;
            s2_ds(:,thickness_ds(2)-n2+1) = a2.*(1-1i*c2);
            s2_ds(:,N2_ds-thickness_ds(2)+n2) = a2.*(1-1i*c2);
            s2_ds(:,N2_ds-thickness_ds(2)+n2) = a2.*(1-1i*c2);

            % sample
            s1_e2 = Device.sampleMap(s1_ds, 'e', 2);
            s1_e3 = Device.sampleMap(s1_ds, 'e', 3);
            s1_h2 = Device.sampleMap(s1_ds, 'h', 2);
            s1_h3 = Device.sampleMap(s1_ds, 'h', 3);
            s2_e1 = Device.sampleMap(s2_ds, 'e', 1);
            s2_e3 = Device.sampleMap(s2_ds, 'e', 3);
            s2_h1 = Device.sampleMap(s2_ds, 'h', 1);
            s2_h3 = Device.sampleMap(s2_ds, 'h', 3);

            % inverse
            S1E2_inv = diag(sparse(1./s1_e2(:)));
            S1E3_inv = diag(sparse(1./s1_e3(:)));
            S1H2_inv = diag(sparse(1./s1_h2(:)));
            S1H3_inv = diag(sparse(1./s1_h3(:)));
            S2E1_inv = diag(sparse(1./s2_e1(:)));
            S2E3_inv = diag(sparse(1./s2_e3(:)));
            S2H1_inv = diag(sparse(1./s2_h1(:)));
            S2H3_inv = diag(sparse(1./s2_h3(:)));
        end

        function ind = setFieldSeq(device_label)
            switch device_label
                case 'xy'
                    ind = [1, 2, 3]; % For 'xy' -> Ex=1, Ey=2, Ez=3
                case 'yz'
                    ind = [2, 3, 1]; % For 'yz' -> Ez=1, Ex=2, Ey=3
                case 'zx'
                    ind = [3, 1, 2]; % For 'xz' -> Ey=1, Ez=2, Ex=3
            end
        end

    end

end