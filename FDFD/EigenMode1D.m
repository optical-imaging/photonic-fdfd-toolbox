classdef (Hidden) EigenMode1D < EigenMode
    %EigenMode1D: 1D FDFD eigenmode solver on an Axis mesh for TE/TM slab cross-sections
    %
    % Key Properties (inherited from FDFD):
    %   lam            - Wavelength (m)
    %   mesh           - Axis mesh used by the eigenmode problem
    %   bc             - BC1D boundary-condition object on the same axis
    %   Ex,Ey,Ez       - Electric-field eigenmodes, size N×Nm
    %   Hx,Hy,Hz       - Magnetic-field eigenmodes, size N×Nm
    %
    % Key Properties (inherited from EigenMode, SetAccess = protected):
    %   modenum        - Requested mode indices (vector)
    %   neff           - Effective indices of solved modes, aligned with modenum
    %
    % Key Properties (EigenMode1D-specific):
    %   pol            - Polarization: 'TE'|'TM' (selects eigen formulation and field components)
    %
    % Dependent Properties:
    %   label          - Axis label from mesh.label when mesh is assigned; [] otherwise
    %   setflag        - True if lam, mesh, bc, and modenum are assigned (inherited dependent)
    %   solflag        - True if all E/H field components are populated (inherited dependent)
    %
    % Key Methods:
    %   EigenMode1D(lam,modenum,pol)     - Construct 1D eigenmode solver with wavelength, modes, and polarization
    %   get.label()                      - Return mesh.label if mesh is assigned
    %
    %   setMesh(axis)                    - Assign Axis mesh to solver (calls FDFD.setMesh)
    %   setBC(bc1d)                      - Assign BC1D and enforce bc label matches solver axis label
    %
    %   solveFDFD(device1d,flipprop)     - Solve eigenmodes on Device1D and populate neff and E/H fields
    %   printInfo()                      - Print wavelength, polarization, modenum, and neff (debug use)
    %
    % Key Methods (Static, Access = private):
    %   assignField(axis_label,...)      - Map local (e1,e2,e3,h1,h2,h3) to (Ex,Ey,Ez,Hx,Hy,Hz) by axis label

    properties
        pol
    end

    properties (Dependent)
        label
    end

    methods
        % constructor
        function obj = EigenMode1D(wavelength, modenum, pol)
            arguments
                wavelength (1,1) double {mustBePositive} % in unit [m]
                modenum (1,:) double {mustBeInteger, mustBePositive}
                pol {mustBeMember(pol, {'TE','TM'})}
            end
            obj@EigenMode(wavelength, modenum);
            obj.pol = pol;
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
        function setMesh(obj, axis)
            arguments
                obj
                axis Axis
            end
            obj.mesh = axis;
        end

        function setBC(obj, bc1d)
            arguments
                obj
                bc1d BC1D
            end
            if ~strcmp(bc1d.label, obj.label)
                dispError('EigenMode1D:BCLabelNotMatch');
            end
            obj.bc = bc1d;
        end

        % launch solver
        function solveFDFD(obj, device, flipprop)
            arguments
                obj
                device Device1D
                flipprop (1,1) logical = false
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
            DE1 = obj.DM1D(obj.mesh, 'E')/k0; % normalized with k0
            DH1 = obj.DM1D(obj.mesh, 'M')/k0;

            % solve eigenmode
            eps = device.eps;
            mu = device.mu;
            switch obj.pol
                case 'TE'
                    % construct coefficient matrix A and solve eigen equation
                    A = DH1*spdiags(1./mu(:,3),0,obj.mesh.n,obj.mesh.n)*DE1 + spdiags(eps(:,2),0,obj.mesh.n,obj.mesh.n);
                    [E2, G2] = eigs(A, spdiags(1./mu(:,1),0,obj.mesh.n,obj.mesh.n), max(obj.modenum), max(eps,[],"all"));
                    neff = sqrt(diag(G2));

                    % reconstruct H fields
                    H1 = -neff.'/eta0.*(spdiags(1./mu(:,1),0,obj.mesh.n,obj.mesh.n)*E2);
                    H3 = 1i/eta0*spdiags(1./mu(:,3),0,obj.mesh.n,obj.mesh.n)*DE1*E2;

                    % sort eigenmodes
                    [~, ind] = sort(real(neff),'descend');
                    e2 = E2(:,ind);
                    h1 = H1(:,ind);
                    h3 = H3(:,ind); % out of phase for Ez and Hy, because wave is propagating towards +x
                    neff = neff(ind);

                    % assign requested modes
                    e2 = e2(:,obj.modenum);
                    h1 = h1(:,obj.modenum);
                    h3 = h3(:,obj.modenum);

                    e1 = zeros(obj.mesh.n, numel(obj.modenum));
                    e3 = zeros(obj.mesh.n, numel(obj.modenum));
                    h2 = zeros(obj.mesh.n, numel(obj.modenum));

                case 'TM'
                    A = DE1*spdiags(1./eps(:,3),0,obj.mesh.n,obj.mesh.n)*DH1+spdiags(mu(:,2),0,obj.mesh.n,obj.mesh.n);
                    [H2, G2] = eigs(A, spdiags(1./eps(:,1),0,obj.mesh.n,obj.mesh.n), max(obj.modenum), max(eps,[],"all"));
                    neff = sqrt(diag(G2));

                    E1 = neff.'*eta0.*(spdiags(1./eps(:,1),0,obj.mesh.n,obj.mesh.n)*H2);
                    E3 = -1i*eta0*spdiags(1./eps(:,3),0,obj.mesh.n,obj.mesh.n)*DH1*H2;

                    [~, ind] = sort(real(neff),'descend');
                    h2 = H2(:,ind);
                    e1 = E1(:,ind);
                    e3 = E3(:,ind); % out of phase for Ez and Hy, because wave is propagating towards +x
                    neff = neff(ind);

                    h2 = h2(:,obj.modenum);
                    e1 = e1(:,obj.modenum);
                    e3 = e3(:,obj.modenum);

                    h1 = zeros(obj.mesh.n, numel(obj.modenum));
                    h3 = zeros(obj.mesh.n, numel(obj.modenum));
                    e2 = zeros(obj.mesh.n, numel(obj.modenum));
            end

            if flipprop
                e2_temp = e2; e3_temp = e3;
                h2_temp = h2; h3_temp = h3;

                e2 = -e3_temp;
                e3 = e2_temp;
                h2 = -h3_temp;
                h3 = h2_temp;
            end

            [E, H] = EigenMode1D.assignField(obj.mesh.label, e1, e2, e3, h1, h2, h3);

            obj.neff = neff(obj.modenum);
            obj.Ex = E.Ex;
            obj.Ey = E.Ey;
            obj.Ez = E.Ez;
            obj.Hx = H.Hx;
            obj.Hy = H.Hy;
            obj.Hz = H.Hz;
        end

        % display
        function printInfo(obj)
            info_name = {
                'Object name';
                'Wavelength';
                'Polarization mode';
                'Mode number';
                'neff'
                };

            lam_str = [num2str(obj.lam), ' m'];
            if isscalar(obj.modenum)
                num_str = num2str(obj.modenum);
            else
                num_str = strjoin(string(obj.modenum), ', ');
            end
            if isscalar(obj.neff)
                neff_str = num2str(obj.neff);
            else
                neff_str = strjoin(string(obj.neff), ', ');
            end

            value = {
                inputname(1);
                lam_str;
                obj.pol;
                num_str;
                neff_str
                };

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-s\n', maxNameLength + 1, info_name{ii}, value{ii});
            end
        end

    end

    methods (Static, Access = private)
        function [E, H] = assignField(axis_label, e1, e2, e3, h1, h2, h3)
            perm = struct('x', [1,2,3], 'y', [3,1,2], 'z', [2,3,1]);
            p = perm.(axis_label);

            E.Ex = eval(sprintf('e%d', p(1))); H.Hx = eval(sprintf('h%d', p(1)));
            E.Ey = eval(sprintf('e%d', p(2))); H.Hy = eval(sprintf('h%d', p(2)));
            E.Ez = eval(sprintf('e%d', p(3))); H.Hz = eval(sprintf('h%d', p(3)));
        end
    end

end