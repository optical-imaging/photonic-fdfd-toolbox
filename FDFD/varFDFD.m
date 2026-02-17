classdef (Hidden) varFDFD < Scattering2D
    %varFDFD: Semi-3D variational FDFD solver using a 3D layer device collapsed to an effective 2D scattering problem
    %
    % Key Properties (inherited from FDFD):
    %   lam            - Wavelength (m)
    %   mesh           - Grid3D mesh for the semi-3D domain (overridden from Scattering2D)
    %   bc             - BC2D boundary-condition container (must be 'xy')
    %   Ex,Ey,Ez       - Reconstructed 3D electric fields, size Nx×Ny×Nz×Ns
    %   Hx,Hy,Hz       - Reconstructed 3D magnetic fields, size Nx×Ny×Nz×Ns
    %
    % Key Properties (inherited from Scattering / Scattering2D):
    %   src            - Source object (Source1D) used for 2D scattering injection
    %   pol            - Total polarization: 'TE'|'TM'
    %
    % Key Properties (varFDFD-specific):
    %   ref            - Reference bundle struct: pin (1×2), refDevice (Device1D), mode (EigenMode1D)
    %   epseff         - Effective 2D index/permittivity map used to build an equivalent Device2D
    %   C              - Integration coefficient bundle struct (C1..C5 and related constants)
    %
    % Dependent Properties:
    %   gridxy         - Grid2D plane extracted from mesh via mesh.getPlane('xy')
    %   label          - Plane label for scattering core (inherited; hard-coded to 'xy' in Scattering2D getter)
    %   setflag        - True if lam, mesh, bc, src, label, ref, epseff, and C are assigned (inherited + extended)
    %   solflag        - True if all E/H field components are populated (inherited dependent)
    %
    % Key Methods (inherited from FDFD/Scattering2D, still used here):
    %   setSrc(src1d)                  - Attach configured Source1D (pol consistency checked in Scattering2D)
    %   solveScattering2DCore(device2d)- Build/split 2D scattering system used as the collapsed core solver
    %
    % Key Methods:
    %   varFDFD(lam,pol)               - Construct variational semi-3D solver (inherits Scattering2D interface)
    %   get.gridxy()                   - Return xy plane Grid2D from Grid3D mesh
    %
    %   setMesh(grid3d)                - Assign Grid3D mesh to solver (overrides Scattering2D signature)
    %   setBC(bc2d)                    - Assign BC2D and enforce bc label is 'xy'
    %
    %   setReference(layerDevice,pin)  - Build reference Device1D + EigenMode1D at pin and compute epseff/C
    %   solveFDFD()                    - Solve collapsed 2D scattering and reconstruct full 3D E/H fields
    %
    % Key Methods (Access = protected):
    %   checkSet()                     - Extend readiness check to require ref, epseff, and C are available
    %
    % Key Methods (Static, Access = private):
    %   refDevice1D(layerDevice,pin)   - Extract 1D z-line Device1D at (x,y)=pin from LayerDevice voxel maps
    %   refEigenMode1D(refDevice,lam,pol) - Solve reference EigenMode1D (collapsed polarization is opposite)
    %   collapseDevice(device,refDev,refMode) - Compute effective 2D map and integration coefficients (C1..C5)
    %   reconstruct3DField([Nx,Ny,Nz],ref_mode,psi) - Reconstruct 3D field as ref_mode(z) * psi(x,y)

    properties
        ref % struct, should include:pin location([1,2]), reference device (Device1D), reference eigenmode (Eigenmode1D)
        epseff % effective 2D map
        C % integration intermidiate constant
    end

    properties (Dependent)
        gridxy
    end

    methods
        % constructor (initialization)
        function obj = varFDFD(wavelength, total_pol)
            arguments
                wavelength double {mustBePositive}
                total_pol {mustBeMember(total_pol, {'TE','TM'})}
            end

            obj@Scattering2D(wavelength, total_pol);
            obj.ref = struct();
            obj.epseff = [];
            obj.C = struct;
        end

        % dependent
        function val = get.gridxy(obj)
            val = obj.mesh.getPlane('xy');
        end

        % set protected properties
        function setMesh(obj, gridxyz)
            arguments
                obj
                gridxyz Grid3D
            end
            obj.mesh = gridxyz;
        end

        function setBC(obj, bc2d)
            arguments
                obj
                bc2d BC2D
            end
            if ~strcmp(bc2d.label, 'xy')
                dispError('varFDFD:BC2DNotXY');
            end
            obj.bc = bc2d;
        end

        % setup reference
        function setReference(obj, device, ref_point)
            arguments
                obj
                device LayerDevice
                ref_point (1,2) double {mustBeReal}
            end

            % find reference eigenmode
            ref_device = obj.refDevice1D(device, ref_point);
            ref_eigenmode = obj.refEigenMode1D(ref_device, obj.lam, obj.pol);

            obj.ref.pin = ref_point;
            obj.ref.refDevice = ref_device;
            obj.ref.mode = ref_eigenmode;

            % C and epseff is depended on ref
            [ind_eff, C_bundle] = obj.collapseDevice(device, ref_device, ref_eigenmode);
            obj.C = C_bundle;
            obj.epseff = ind_eff;
        end

        % varFDFD propagation
        function solveFDFD(obj)
            % simulation device is epseff
            % pre-check
            if ~obj.setflag
                dispError('varFDFD:setvarFDFDFirst');
            end

            % constant
            k0 = 2*pi/obj.lam; % wave number
            beta = k0*obj.src.neff;
            Nx = obj.mesh.axisx.n;
            Ny = obj.mesh.axisy.n;
            Nz = obj.mesh.axisz.n;
            dx = obj.mesh.axisx.d;
            dy = obj.mesh.axisy.d;

            % create effective device2D
            Alpha = Device2D(1);
            eps_eff = obj.epseff;
            mu_eff = ones(size(eps_eff));
            Alpha.setBitmap(eps_eff,mu_eff,[dx dy],obj.gridxy.label);
            Alpha.shiftMesh('x', 'min', obj.gridxy.axis1.v(1)); % align mesh
            Alpha.shiftMesh('y', 'min', obj.gridxy.axis2.v(1));

            % solve 2D scattering problem
            [psi_2d, D_pack, scpml_pack, ~, ~] = obj.solveScattering2DCore(Alpha);

            % reference integration coefficients
            ref_mode = obj.ref.mode;
            C1 = obj.C.C1; C2 = obj.C.C2;
            C3 = obj.C.C3; C4 = obj.C.C4;

            % reconstruction
            EPS_EFF_inv = spdiags(1./obj.epseff(:), 0, Nx*Ny, Nx*Ny);
            C3_diag_inv = spdiags(1./C3(:),0,Nx*Ny,Nx*Ny);
            C4_diag = spdiags(C4(:),0,Nx*Ny,Nx*Ny);
            n_fields = size(psi_2d,2);

            Ex = zeros(Nx, Ny, Nz, n_fields);
            Ey = zeros(Nx, Ny, Nz, n_fields);
            Ez = zeros(Nx, Ny, Nz, n_fields);
            Hx = zeros(Nx, Ny, Nz, n_fields);
            Hy = zeros(Nx, Ny, Nz, n_fields);
            Hz = zeros(Nx, Ny, Nz, n_fields);

            for kk = 1:n_fields
                switch obj.pol
                    case 'TE'
                        % x-y psi pattern
                        psi_Ez = psi_2d(:,kk);
                        psi_Hx = (1i*beta*C4_diag)/(k0*C2)*EPS_EFF_inv*scpml_pack.S2H1_inv*D_pack.DE2*psi_Ez; % DE contains 1/k0
                        psi_Hy = -(1i*beta*C4_diag)/(k0*C2)*EPS_EFF_inv*scpml_pack.S1H2_inv*D_pack.DE1*psi_Ez;
                        psi_Ex = -C1*C3_diag_inv*psi_Hy;
                        psi_Ey = C1*C3_diag_inv*psi_Hx;

                        % reference mode
                        hx = ref_mode.Hx; hy = hx;
                        ey = ref_mode.Ey; ex = ey;
                        ez = ref_mode.Ez;

                        % reconstruct
                        Ez(:,:,:,kk) = varFDFD.reconstruct3DField([Nx,Ny,Nz],ez,psi_Ez);
                        Hz(:,:,:,kk) = zeros(size(Ez));

                    case 'TM'
                        % x-y psi pattern
                        psi_Hz = psi_2d(:,kk);
                        psi_Ex = (1i*beta*C4_diag)/(k0*C2)*EPS_EFF_inv*scpml_pack.S2E1_inv*D_pack.DH2*psi_Hz; % DH contains 1/k0
                        psi_Ey = -(1i*beta*C4_diag)/(k0*C2)*EPS_EFF_inv*scpml_pack.S1E2_inv*D_pack.DH1*psi_Hz;
                        psi_Hx = -C1*C3_diag_inv*psi_Ey;
                        psi_Hy = C1*C3_diag_inv*psi_Ex;

                        % reference mode
                        ex = ref_mode.Ex; ey = ex;
                        hy = ref_mode.Hy; hx = hy;
                        hz = ref_mode.Hz;

                        % reconstruct
                        Hz(:,:,:,kk) = varFDFD.reconstruct3DField([Nx,Ny,Nz],hz,psi_Hz);
                        Ez(:,:,:,kk) = zeros(size(Ez));
                end

                Ex(:,:,:,kk) = varFDFD.reconstruct3DField([Nx,Ny,Nz],ex,psi_Ex);
                Ey(:,:,:,kk) = varFDFD.reconstruct3DField([Nx,Ny,Nz],ey,psi_Ey);
                Hx(:,:,:,kk) = varFDFD.reconstruct3DField([Nx,Ny,Nz],hx,psi_Hx);
                Hy(:,:,:,kk) = varFDFD.reconstruct3DField([Nx,Ny,Nz],hy,psi_Hy);
            end

            obj.Ex = Ex;
            obj.Ey = Ey;
            obj.Ez = Ez;
            obj.Hx = Hx;
            obj.Hy = Hy;
            obj.Hz = Hz;
        end


    end

    methods (Access = protected)
        function temp_status = checkSet(obj)
            temp_status = obj.checkSet@Scattering2D;
            temp_status = temp_status && ~isempty(obj.ref)...
                && ~isempty(obj.epseff)...
                && ~isempty(obj.C);
        end
    end

    methods(Static, Access = private)
        function ref_device = refDevice1D(device, ref_point)
            ref_device = Device1D(1); % empty device

            % find ref pin locatopm
            [~, x_idx] = min(abs(device.mesh.axisx.v - ref_point(1)));
            [~, y_idx] = min(abs(device.mesh.axisy.v - ref_point(2)));

            eps_ref = squeeze(device.eps(x_idx,y_idx,:,3));
            mu_ref = squeeze(device.mu(x_idx,y_idx,:,3));
            zaxis = device.mesh.axisz;

            % directly assign bitmap
            ref_device.setBitmap(eps_ref,mu_ref,zaxis.d,'z');
            ref_device.mesh.shiftGrid('min', zaxis.v(1));
        end

        function ref_eigenmode = refEigenMode1D(ref_device,wavelength,total_pol)
            ref_axis = ref_device.mesh;
            bc_ref = BC1D('z', 'Perfect'); % Perfect boundary condition along z (to be collapsed)

            % create and set up EigenMode1D object
            switch total_pol
                case 'TE'
                    ref_eigenmode = EigenMode1D(wavelength, 1, 'TM'); % entire TE, then the collapsed eigenmode is TM
                case 'TM'
                    ref_eigenmode = EigenMode1D(wavelength, 1, 'TE'); % entire TM, then the collapsed eigenmode is TE
            end
            ref_eigenmode.setMesh(ref_axis);
            ref_eigenmode.setBC(bc_ref);

            % solve 1D reference eigenmode
            ref_eigenmode.solveFDFD(ref_device, false);
        end

        function [ind_eff, C] = collapseDevice(device, ref_device, ref_eigenmode)
            % no 2x grid at this stage
            % prepare
            Nx = device.mesh.axisx.n;
            Ny = device.mesh.axisy.n;
            Nz = device.mesh.axisz.n;

            eps_xyz = device.eps(:,:,:,3); eps_ref = ref_device.eps(:,1);
            mu_xyz = device.mu(:,:,:,3); mu_ref = ref_device.mu(:,1);

            ref_mesh = ref_device.mesh;

            k0 = 2*pi/ref_eigenmode.lam;

            % varFDFD integration coefficients
            switch ref_eigenmode.pol % the entire polarization should be opposite
                case 'TM' % for total TE
                    dhz = FDFD.DM1D(ref_mesh,'M');
                    hx_mode = ref_eigenmode.Hx;

                    hx_pad = repmat(reshape(hx_mode, [1 1 Nz]),[Nx Ny 1]);
                    diff_hx_pad = repmat(reshape(dhz*hx_mode, [1 1 Nz]),[Nx Ny 1]);
                    eps_ref_pad = repmat(reshape(eps_ref, [1 1 Nz]),[Nx Ny 1]);

                    C1 = sum((dhz*hx_mode).^2./eps_ref);
                    C2 = sum(hx_mode.^2./eps_ref);
                    C3 = sum( eps_xyz./eps_ref_pad.^2.*diff_hx_pad.^2, 3);
                    C4 = sum( eps_xyz./eps_ref_pad.^2.*hx_pad.^2, 3);
                    C5 = sum( mu_xyz.*hx_pad.^2, 3);

                case 'TE'

                    dez = FDFD.DM1D(ref_mesh,'E');
                    ex_mode = ref_eigenmode.Ex;

                    ex_pad = repmat(reshape(ex_mode, [1 1 Nz]),[Nx Ny 1]);
                    diff_ex_pad = repmat(reshape(dez*ex_mode, [1 1 Nz]),[Nx Ny 1]);
                    mu_ref_pad = repmat(reshape(mu_ref, [1 1 Nz]),[Nx Ny 1]);

                    C1 = sum((dez*ex_mode).^2./mu_ref);
                    C2 = sum(ex_mode.^2./mu_ref);
                    C3 = sum( mu_xyz./mu_ref_pad.^2.*diff_ex_pad.^2, 3);
                    C4 = sum( mu_xyz./mu_ref_pad.^2.*ex_pad.^2, 3);
                    C5 = sum( eps_xyz.*ex_pad.^2, 3);

            end

            ind_eff = C4/C2^2.*(C5-C1^2./(k0^2*C3));
            eff_th = 1;
            ind_eff(ind_eff<eff_th) = eff_th;

            C = struct();
            C.C1 = C1;
            C.C2 = C2;
            C.C3 = C3;
            C.C4 = C4;
            C.C5 = C5;
        end

        function full_field = reconstruct3DField(N, ref_mode, psi)
            Nx = N(1); Ny = N(2); Nz = N(3);

            psi = reshape(psi,[Nx Ny 1]);
            ref_mode = reshape(ref_mode,[1 1 Nz]);
            ref_mode_rep = repmat(ref_mode,[Nx Ny 1]);

            full_field = ref_mode_rep.*psi;
        end

    end



end