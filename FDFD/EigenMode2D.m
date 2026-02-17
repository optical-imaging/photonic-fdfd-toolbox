classdef (Hidden) EigenMode2D < EigenMode
    %EigenMode2D: 2D FDFD eigenmode solver on a Grid2D mesh for waveguide cross-sections
    %
    % Key Properties (inherited from FDFD):
    %   lam            - Wavelength (m)
    %   mesh           - Grid2D mesh used by the eigenmode problem
    %   bc             - BC2D boundary-condition object on the same plane
    %   Ex,Ey,Ez       - Electric-field eigenmodes, size N1×N2×Nm
    %   Hx,Hy,Hz       - Magnetic-field eigenmodes, size N1×N2×Nm
    %
    % Key Properties (inherited from EigenMode, SetAccess = protected):
    %   modenum        - Requested mode indices (vector)
    %   neff           - Effective indices of solved modes, aligned with modenum
    %
    % Dependent Properties:
    %   label          - Plane label from mesh.label when mesh is assigned; [] otherwise
    %   setflag        - True if lam, mesh, bc, and modenum are assigned (inherited dependent)
    %   solflag        - True if all E/H field components are populated (inherited dependent)
    %
    % Key Methods:
    %   EigenMode2D(lam,modenum)         - Construct 2D eigenmode solver with wavelength and mode indices
    %   get.label()                       - Return mesh.label if mesh is assigned
    %
    %   setMesh(grid2d)                   - Assign Grid2D mesh to solver (calls FDFD.setMesh)
    %   setBC(bc2d)                       - Assign BC2D and enforce bc label matches solver plane label
    %
    %   solveFDFD(device2d)               - Solve 2D eigenmodes on Device2D and populate neff and E/H fields
    %   printInfo()                       - Print wavelength, modenum, and neff (debug use)

    properties (Dependent)
        label
    end

    methods
        % constructor
        function obj = EigenMode2D(wavelength, modenum)
            arguments
                wavelength (1,1) double {mustBePositive} % in unit [m]
                modenum (1,:) double {mustBeInteger, mustBePositive}
            end
            obj@EigenMode(wavelength, modenum);
        end

        % dependent
        function val = get.label(obj)
            if ~isempty(obj.mesh)
                val = obj.mesh.label;
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
                dispError('EigenMode1D:BCLabelNotMatch');
            end
            obj.bc = bc2d;
        end

        % launch solver
        function solveFDFD(obj, device)
            arguments
                obj
                device Device2D
            end

            % pre-check
            if ~obj.setflag
                dispError('EigenMode:SetUpEigenSolver');
            end
            if ~device.meshflag
                dispError('Device2D:MeshDevice');
            end

            % Constant
            eta0 = Constant("eta0").v;
            k0 = 2*pi/obj.lam;

            % derivative matrices
            [DE1, DE2] = obj.DM2D(obj.mesh, 'E');
            [DH1, DH2] = obj.DM2D(obj.mesh, 'M');
            DE1 = DE1/k0;
            DE2 = DE2/k0;
            DH1 = DH1/k0;
            DH2 = DH2/k0;

            % device meshed array
            eps = device.eps;
            mu = device.mu;
            N1 = obj.mesh.axis1.n;
            N2 = obj.mesh.axis2.n;

            % diagonize eps and mu map
            eps1 = eps(:,:,1); eps1_inv = spdiags(1./eps1(:),0,N1*N2,N1*N2); eps1 = spdiags(eps1(:),0,N1*N2,N1*N2); % eps_xx (conventional)
            eps2 = eps(:,:,2); eps2_inv = spdiags(1./eps2(:),0,N1*N2,N1*N2); eps2 = spdiags(eps2(:),0,N1*N2,N1*N2); % eps_yy
            eps3 = eps(:,:,3); eps3_inv = spdiags(1./eps3(:),0,N1*N2,N1*N2); eps3 = spdiags(eps3(:),0,N1*N2,N1*N2); % eps_zz
            mu1 = mu(:,:,1); mu1_inv = spdiags(1./mu1(:),0,N1*N2,N1*N2); mu1 = spdiags(mu1(:),0,N1*N2,N1*N2); % mu_xx
            mu2 = mu(:,:,2); mu2_inv = spdiags(1./mu2(:),0,N1*N2,N1*N2); mu2 = spdiags(mu2(:),0,N1*N2,N1*N2); % mu_yy
            mu3 = mu(:,:,3); mu3_inv = spdiags(1./mu3(:),0,N1*N2,N1*N2); mu3 = spdiags(mu3(:),0,N1*N2,N1*N2); % mu_zz

            % PQ matrix
            P = [ DE1*eps3_inv*DH2, -(DE1*eps3_inv*DH1+mu2);...
                DE2*eps3_inv*DH2+mu1, -DE2*eps3_inv*DH1 ];
            Q = [ DH1*mu3_inv*DE2, -(DH1*mu3_inv*DE1+eps2);...
                DH2*mu3_inv*DE2+eps1, -DH2*mu3_inv*DE1];

            % eigenmode problem
            eps_core = max(eps,[],"all","ComparisonMethod","real");
            [E,G2] = eigs(P*Q, max(obj.modenum), -eps_core);
            neff = sqrt(-diag(G2));

            % calculate H field
            H = zeros(size(E));
            for ii = 1:max(obj.modenum)
                H(:,ii) = 1/(eta0*neff(ii))*Q*E(:,ii);
            end

            % extract transverse parts
            e1 = E(1:N1*N2,:);
            e2 = E(N1*N2+1:end,:);
            h1 = H(1:N1*N2,:);
            h2 = H(N1*N2+1:end,:);

            % calculate longitude components
            e3 = -1i*eta0*eps3_inv*(DH1*h2 - DH2*h1);
            h3 = 1i/eta0*mu3_inv*(DE1*e2 - DE2*e1);

            % reshape
            E1 = zeros(N1, N2, numel(obj.modenum));
            E2 = zeros(N1, N2, numel(obj.modenum));
            E3 = zeros(N1, N2, numel(obj.modenum));
            H1 = zeros(N1, N2, numel(obj.modenum));
            H2 = zeros(N1, N2, numel(obj.modenum));
            H3 = zeros(N1, N2, numel(obj.modenum));
            for ii = 1:numel(obj.modenum)
                E1(:,:,ii) = reshape(e1(:,obj.modenum(ii)),[N1, N2]);
                E2(:,:,ii) = reshape(e2(:,obj.modenum(ii)),[N1, N2]);
                E3(:,:,ii) = reshape(e3(:,obj.modenum(ii)),[N1, N2]);
                H1(:,:,ii) = reshape(h1(:,obj.modenum(ii)),[N1, N2]);
                H2(:,:,ii) = reshape(h2(:,obj.modenum(ii)),[N1, N2]);
                H3(:,:,ii) = reshape(h3(:,obj.modenum(ii)),[N1, N2]);
            end

            % assaign field components
            Ex = zeros(size(E1));
            Ey = zeros(size(E1));
            Ez = zeros(size(E1));
            Hx = zeros(size(E1));
            Hy = zeros(size(E1));
            Hz = zeros(size(E1));

            eval(sprintf('E%s = E1;', obj.label(1)));
            eval(sprintf('E%s = E2;', obj.label(2)));
            eval(sprintf('H%s = H1;', obj.label(1)));
            eval(sprintf('H%s = H2;', obj.label(2)));

            long_axis = setdiff('xyz', obj.label);
            eval(sprintf('E%s = E3;', long_axis));
            eval(sprintf('H%s = H3;', long_axis));

            obj.neff = neff(obj.modenum);
            obj.Ex = Ex;
            obj.Ey = Ey;
            obj.Ez = Ez;
            obj.Hx = Hx;
            obj.Hy = Hy;
            obj.Hz = Hz;
        end

        % Display
        function printInfo(obj)
            info_name = {'Object Name';
                'Wavelength';
                'Mode';
                'neff'};
            value = {inputname(1);
                [num2str(obj.lam), ' m'];
                obj.modenum;
                obj.neff};

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength+1, info_name{ii}, value{ii});
            end
        end

    end

end