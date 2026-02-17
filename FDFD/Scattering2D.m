classdef (Hidden) Scattering2D < Scattering
    %Scattering2D: 2D FDFD scattering solver with 1D port source injection and PML/PBC support
    %
    % Key Properties (inherited from FDFD):
    %   lam            - Wavelength (m)
    %   mesh           - Grid2D mesh for the scattering domain (or Grid3D when used by varFDFD path)
    %   bc             - BC2D boundary-condition container for the 2D domain
    %   Ex,Ey,Ez       - Solved electric-field components, size N1×N2×Ns
    %   Hx,Hy,Hz       - Solved magnetic-field components, size N1×N2×Ns
    %
    % Key Properties (inherited from Scattering):
    %   src            - Source object providing incident fields (must be Source1D here)
    %
    % Key Properties (Scattering2D-specific):
    %   pol            - Polarization: 'TE'|'TM' (selects scalar solve variable and reconstruction)
    %
    % Constant Properties (SCPML parameters):
    %   scamp          - PML amplitude scaling factor
    %   conduct        - PML conductivity scaling factor
    %   power          - PML polynomial order
    %
    % Dependent Properties:
    %   label          - Plane label derived from mesh (special-cased for varFDFD: returns 'xy')
    %   setflag        - True if lam, mesh, bc, src, and label are assigned (inherited + extended)
    %   solflag        - True if all E/H field components are populated (inherited dependent)
    %
    % Key Methods (inherited from FDFD/Scattering, still used here):
    %   Scattering(lam)              - Construct base scattering solver
    %   setMesh(mesh)                - Assign mesh (overridden here for Grid2D signature)
    %   setBC(bc)                    - Assign boundary conditions (overridden here for BC2D signature)
    %   setSrc(src)                  - Attach configured Source (overridden here for Source1D + pol check)
    %   get.solflag()                - Return solved flag
    %
    % Key Methods:
    %   Scattering2D(lam,pol)        - Construct 2D scattering solver with wavelength and polarization
    %   get.label()                  - Return mesh plane label (or 'xy' for varFDFD)
    %   setMesh(grid2d)              - Assign Grid2D mesh to solver
    %   setBC(bc2d)                  - Assign BC2D and enforce bc label matches solver label
    %   setSrc(src1d)                - Attach Source1D and enforce polarization match
    %   solveFDFD(device2d)          - Solve scattering problem and reconstruct all E/H components
    %
    % Key Methods (Access = protected):
    %   checkSet()                   - Extend readiness check to require label is available
    %   solveScattering2DCore(device)- Build system matrix, source RHS, and solve for scalar field vector(s)
    %
    % Key Methods (Static, Access = protected):
    %   rotateProfile(...)           - Rotate/interpolate source transverse profile for oblique incidence
    %   rotatePhase(...)             - Generate rotated propagation phase pattern on 2x grid
    %   SCPML(N1ds,N2ds,PML)         - Build stretched-coordinate PML scaling matrices (inverse diagonals)
    %
    % Key Methods (Static, Access = private):
    %   setFieldSeq(label)           - Map (e1,e2,e3)/(h1,h2,h3) into (Ex,Ey,Ez)/(Hx,Hy,Hz) ordering by plane

    properties
        pol
    end

    properties (Constant) % SCPML for 2D
        scamp = 3;
        conduct = 150;
        power = 3;
    end

    properties (Dependent)
        label
    end

    methods
        % constructor
        function obj = Scattering2D(wavelength, pol)
            arguments
                wavelength (1,1) double {mustBePositive} % in unit [m]
                pol {mustBeMember(pol, {'TE','TM'})}
            end
            obj@Scattering(wavelength);
            obj.pol = pol;
        end

        % dependent
        function val = get.label(obj)
            if ~isempty(obj.mesh)
                switch class(obj)
                    case 'Scattering2D'
                        val = obj.mesh.label;
                    case 'varFDFD'
                        val = 'xy';
                end
            else
                val = [];
            end
        end

        % set protected properties
        function setMesh(obj, mesh2d)
            arguments
                obj
                mesh2d Grid2D
            end
            obj.mesh = mesh2d;
        end

        function setBC(obj, bc2d)
            arguments
                obj
                bc2d BC2D
            end
            if ~strcmp(bc2d.label, obj.label)
                dispError('Scattering2D:BC2DLabelNotMatch');
            end
            obj.bc = bc2d;
        end

        % set source
        function setSrc(obj,src)
            arguments
                obj
                src Source1D
            end
            if ~strcmp(obj.pol, src.pol)
                dispError();
            end
            obj.setSrc@Scattering(src);
        end

        % launch solver
        function solveFDFD(obj, device)
            arguments
                obj
                device Device2D
            end

            % pre-check
            if ~obj.setflag
                dispError('Scattering2D:ScatteringNotSetUp');
            end

            [f, D_pack, scpml_pack, eps_inv, mu_inv] = obj.solveScattering2DCore(device);

            % constant
            eta0 = Constant('eta0').v; % vacuum impedance
            n_fields = size(f,2);

            % reconstruct and assign all fields
            e1 = zeros(obj.mesh.axis1.n*obj.mesh.axis2.n, n_fields);
            e2 = zeros(obj.mesh.axis1.n*obj.mesh.axis2.n, n_fields);
            e3 = zeros(obj.mesh.axis1.n*obj.mesh.axis2.n, n_fields);
            h1 = zeros(obj.mesh.axis1.n*obj.mesh.axis2.n, n_fields);
            h2 = zeros(obj.mesh.axis1.n*obj.mesh.axis2.n, n_fields);
            h3 = zeros(obj.mesh.axis1.n*obj.mesh.axis2.n, n_fields);

            switch obj.pol
                case 'TE'
                    e3 = f;
                    h1 = mu_inv.mu1_inv*scpml_pack.S2H1_inv*D_pack.DE2*e3;
                    h2 = -mu_inv.mu2_inv*scpml_pack.S1H2_inv*D_pack.DE1*e3;
                case 'TM'
                    h3 = f;
                    e1 = eps_inv.eps1_inv*scpml_pack.S2E1_inv*D_pack.DH2*h3;
                    e2 = -eps_inv.eps2_inv*scpml_pack.S1E2_inv*D_pack.DH1*h3;
            end

            % assign fields
            E_field_list = {'Ex', 'Ey', 'Ez'};
            H_field_list = {'Hx', 'Hy', 'Hz'};
            field_seq = Scattering2D.setFieldSeq(obj.mesh.label);
            obj.(E_field_list{field_seq(1)}) = reshape(e1,...
                [obj.mesh.axis1.n, obj.mesh.axis2.n, n_fields]);
            obj.(H_field_list{field_seq(1)}) = reshape(h1,...
                [obj.mesh.axis1.n, obj.mesh.axis2.n, n_fields])/(-1i*eta0);
            obj.(E_field_list{field_seq(2)}) = reshape(e2,...
                [obj.mesh.axis1.n, obj.mesh.axis2.n, n_fields]);
            obj.(H_field_list{field_seq(2)}) = reshape(h2,...
                [obj.mesh.axis1.n, obj.mesh.axis2.n, n_fields])/(-1i*eta0);
            obj.(E_field_list{field_seq(3)}) = reshape(e3,...
                [obj.mesh.axis1.n, obj.mesh.axis2.n, n_fields]);
            obj.(H_field_list{field_seq(3)}) = reshape(h3,...
                [obj.mesh.axis1.n, obj.mesh.axis2.n, n_fields])/(-1i*eta0);
        end

    end

    methods (Access = protected)
        function temp_status = checkSet(obj)
            temp_status = obj.checkSet@Scattering;
            temp_status = temp_status && ~isempty(obj.label);
        end

        function [field_vec, D_pack, scpml_pack, eps_inv, mu_inv] = solveScattering2DCore(obj, device)
            arguments
                obj
                device Device2D
            end

            % Note: obj.mesh is Grid2D for Scattering2D but Grid3D for
            % varFDFD. To be consistant, use device's mesh below
            % For varFDFD, device is effective device made of epseff
            if ~device.meshflag
                dispError('Device2D:MeshDevice');
            end
            mesh2d = device.mesh;

            k0 = 2*pi/obj.lam; % wave number

            % source injection direction
            ind_axis_port = find(mesh2d.label == obj.src.port.dir(2));
            switch ind_axis_port
                case 1
                    axis_source = mesh2d.axis2;
                    axis_prop = mesh2d.axis1;
                    mesh_tp_flag = false;
                case 2
                    axis_source = mesh2d.axis1;
                    axis_prop = mesh2d.axis2;
                    mesh_tp_flag = true;
            end

            % prepare source
            extended_src_field = obj.src.extendField(axis_source);
            field_pol = setdiff('xyz', mesh2d.label, 'stable');
            switch obj.pol
                case 'TE'
                    source_field = extended_src_field.(['E',field_pol]);
                case 'TM'
                    source_field = extended_src_field.(['H',field_pol]);
            end

            % 2x grid
            axis1_ds = mesh2d.axis1.doublesampleAxis;
            axis2_ds = mesh2d.axis2.doublesampleAxis;
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
            [DE1, DE2] = obj.DM2D(mesh2d, 'E'); % PML boundary condition as default
            [DH1, DH2] = obj.DM2D(mesh2d, 'M');

            if any(strcmp({obj.bc.bc1.type, obj.bc.bc2.type}, 'Periodic')) % check periodic condition

                switch mesh_tp_flag
                    case false % incidence towards axis 1
                        k_inc = obj.src.neff*k0*[cos(obj.src.theta), sin(obj.src.theta)];
                    case true % incidence towards axis 2
                        k_inc = obj.src.neff*k0*[sin(obj.src.theta), cos(obj.src.theta)];
                end

                pbc_label = strcmp(obj.bc.bc1.type, 'Periodic') + 2*strcmp(obj.bc.bc2.type, 'Periodic');
                switch pbc_label
                    case 1 % axis 1 boundary
                        [DE1, ~] = obj.DM2D_PBC(mesh2d, 'E', k_inc);
                        [DH1, ~] = obj.DM2D_PBC(mesh2d, 'M', k_inc);
                    case 2 % axis 2 boundary
                        [~, DE2] = obj.DM2D_PBC(mesh2d, 'E', k_inc);
                        [~, DH2] = obj.DM2D_PBC(mesh2d, 'M', k_inc);
                    case 3 % both
                        [DE1, DE2] = obj.DM2D_PBC(mesh2d, 'E', k_inc);
                        [DH1, DH2] = obj.DM2D_PBC(mesh2d, 'M', k_inc);
                end
            end

            DE1 = DE1/k0;
            DE2 = DE2/k0;
            DH1 = DH1/k0;
            DH2 = DH2/k0;

            % SCPML
            PML.scamp = obj.scamp;
            PML.conduct = obj.conduct;
            PML.power = obj.power;
            PML.thickness = [obj.bc.bc1.pmlnum, obj.bc.bc2.pmlnum];
            [S1E2_inv, S1E3_inv, S1H2_inv, S1H3_inv,...
                S2E1_inv, S2E3_inv, S2H1_inv, S2H3_inv] ...
                = Scattering2D.SCPML(2*mesh2d.axis1.n, 2*mesh2d.axis2.n, PML);

            % A matrix
            switch obj.pol
                case 'TE'
                    mu1_inv = sparse(1./device.mu(:,:,1)); mu1_inv = diag(mu1_inv(:));
                    mu2_inv = sparse(1./device.mu(:,:,2)); mu2_inv = diag(mu2_inv(:));
                    eps3 = device.eps(:,:,3); eps3 = diag(sparse(eps3(:)));
                    A = S1E3_inv*DH1*mu2_inv*S1H2_inv*DE1 + S2E3_inv*DH2*mu1_inv*S2H1_inv*DE2 + eps3;

                    mu_inv = struct();
                    mu_inv.mu1_inv = mu1_inv;
                    mu_inv.mu2_inv = mu2_inv;
                    eps_inv = struct();
                case 'TM'
                    eps1_inv = sparse(1./device.eps(:,:,1)); eps1_inv = diag(eps1_inv(:));
                    eps2_inv = sparse(1./device.eps(:,:,2)); eps2_inv = diag(eps2_inv(:));
                    mu3 = device.mu(:,:,3); mu3 = diag(sparse(mu3(:)));
                    A = S1H3_inv*DE1*eps2_inv*S1E2_inv*DH1 + S2H3_inv*DE2*eps1_inv*S2E1_inv*DH2 + mu3;

                    eps_inv = struct();
                    eps_inv.eps1_inv = eps1_inv;
                    eps_inv.eps2_inv = eps2_inv;
                    mu_inv = struct();
            end

            % Q masking matrix
            Q = true(mesh2d.axis1.n, mesh2d.axis2.n);
            switch mesh_tp_flag
                case false
                    [~, ind_p] = min(abs(axis_prop.v-obj.src.port.p1(1)));
                    Q(ind_p:end,:) = false;
                case true
                    [~, ind_p] = min(abs(axis_prop.v-obj.src.port.p1(2)));
                    Q(:,ind_p:end) = false;
            end

            switch obj.src.dir(1)
                case '+'
                    Q = diag(sparse(Q(:)));
                case '-'
                    Q = diag(sparse(~Q(:)));
            end

            % incorporate source field
            n_source = numel(obj.src.neff); % source defined on the same port
            b = zeros(mesh2d.axis1.n*mesh2d.axis2.n, n_source);
            for ii = 1:n_source
                fsrc_profile = Scattering2D.rotateProfile(obj.src.port, mesh_tp_flag, source_field(:,ii), axis_source, axis_prop, obj.src.theta);
                if mesh_tp_flag
                    fsrc_profile = fsrc_profile.';
                end

                fsrc_phase_ds = Scattering2D.rotatePhase(obj.src.neff(ii)*k0,...
                    V_source, V_prop, axis_prop.v(ind_p), obj.src.theta, obj.src.phi);
                switch obj.src.pol
                    case 'TE'
                        fsrc = fsrc_profile.*fsrc_phase_ds(1:2:end, 1:2:end);
                    case 'TM'
                        fsrc = fsrc_profile.*fsrc_phase_ds(2:2:end, 2:2:end);
                end

                b(:,ii) = (Q*A-A*Q)*fsrc(:);
            end

            % solve field
            field_vec = zeros(mesh2d.axis1.n*mesh2d.axis2.n, n_source);
            for kk = 1:n_source
                b_sub = b(:,kk);
                field_vec(:,kk) = A\b_sub;
            end

            % pack useful matrices for reconstruction
            D_pack = struct();
            D_pack.DE1 = DE1;
            D_pack.DE2 = DE2;
            D_pack.DH1 = DH1;
            D_pack.DH2 = DH2;

            scpml_pack = struct();
            scpml_pack.S1E2_inv = S1E2_inv;
            scpml_pack.S1E3_inv = S1E3_inv;
            scpml_pack.S1H2_inv = S1H2_inv;
            scpml_pack.S1H3_inv = S1H3_inv;
            scpml_pack.S2E1_inv = S2E1_inv;
            scpml_pack.S2E3_inv = S2E3_inv;
            scpml_pack.S2H1_inv = S2H1_inv;
            scpml_pack.S2H3_inv = S2H3_inv;
        end
    end

    methods (Static, Access = protected)
        % source field profile rotation
        function profile_rot = rotateProfile(port, mesh_tp_flag, source_field, axis_source, axis_prop, theta)
            arguments
                port Port1D
                mesh_tp_flag
                source_field (:,1)
                axis_source
                axis_prop
                theta
            end

            % plane wave case
            if all(source_field(:) == source_field(1))
                profile_rot = repmat(source_field, [1, axis_prop.n]);
                return
            end

            % rotate grid
            theta = theta/180*pi;
            [PROP, SOURCE] = ndgrid(axis_prop.v, axis_source.v);

            switch mesh_tp_flag
                case false
                    prop_c = (port.p1(1)+port.p2(1))/2;
                    source_c = (port.p1(2)+port.p2(2))/2;
                case true
                    prop_c = (port.p1(2)+port.p2(2))/2;
                    source_c = (port.p1(1)+port.p2(1))/2;
            end

            PROPc = PROP-prop_c;
            SOURCEc = SOURCE-source_c;

            [TH, R] = cart2pol(PROPc, SOURCEc);
            [~, SOURCEc] = pol2cart(TH-theta, R);

            SOURCE = SOURCEc + source_c;

            % interpolate source
            profile_interp = griddedInterpolant(axis_source.v, source_field,...
                'linear', 'nearest');
            profile_rot = profile_interp(SOURCE);
            profile_rot(isnan(profile_rot)) = 0;
        end

        % source injection propagation phase pattern rotation
        function phase2_rot = rotatePhase(beta, AXIS_ds_source, AXIS_ds_prop, port_position, theta, phi)
            % rotate grid
            theta = theta/180*pi;
            PROP = AXIS_ds_prop-port_position;
            SOURCE = AXIS_ds_source;
            [TH, R] = cart2pol(PROP, SOURCE);
            [PROP, ~] = pol2cart(TH-theta, R);

            % assign new phase pattern
            phase2_rot = exp(-1i*beta*PROP+1i*phi);
        end

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
    end

    methods (Static, Access = private)
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