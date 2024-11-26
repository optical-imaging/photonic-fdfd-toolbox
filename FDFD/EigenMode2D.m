classdef (Hidden) EigenMode2D < FDFD2D
    % EigenMode2D: 2D eigenmode solver for electromagnetic analysis
    %
    % Description:
    %   The `EigenMode2D` class is used for computing the eigenmodes of two-dimensional electromagnetic structures.
    %   It provides capabilities to solve for the effective refractive index and the field distributions (electric and magnetic fields).
    %   This class inherits from `FDFD2D` and utilizes finite-difference frequency-domain (FDFD) techniques for calculations.
    %
    % Properties:
    %   mesh - The defined mesh grid for the 2D domain.
    %   lam - The wavelength of the electromagnetic wave in the simulation.
    %   solflag - A flag indicating if the solution has been computed.
    %   Ex, Ey, Ez - Electric field components along x, y, and z dimensions.
    %   Hx, Hy, Hz - Magnetic field components along x, y, and z dimensions.
    %   num - Mode number to solve for.
    %   neff - Effective refractive index of the computed mode.
    %
    % Methods:
    %   EigenMode2D(eps, mu, mesh, plane_label, wavelength, num_mode) - Constructor to create a 2D eigenmode solver object.
    %   printInfo() - Displays information about the computed mode.
    %   dispImg(plot_part) - Displays the eigenmode field for visualization.
    %
    % Example:
    %   % Define material properties, mesh grid, and wavelength
    %   eps = rand(100, 100, 3); % Example relative permittivity map
    %   mu = ones(100, 100, 3); % Relative permeability map
    %   mesh = Grid2D(Axis('x', 0, 1, 100), Axis('y', 0, 1, 100));
    %   plane_label = 'xy';
    %   wavelength = 1.55e-6;
    %
    %   % Create an EigenMode2D solver object and compute the eigenmode
    %   emode = EigenMode2D(eps, mu, mesh, plane_label, wavelength, 1);
    %   emode.printInfo();
    %
    % Notes:
    %   - This class is marked as "Hidden", which means it is not intended to be directly accessed by users.
    %   - Ensure that the `Grid2D` object is properly defined before using it.
    %
    % See Also:
    %   FDFD, FDFD2D, Grid, Axis, Grid2D, Device2D

    properties (SetAccess = private)
        num
        neff
    end

    methods

        function obj = EigenMode2D(eps, mu, mesh, plane_label, wavelength, num_mode)
            % EigenMode2D: Constructor to create a 2D eigenmode solver object
            %
            %   Syntax:
            %     obj = EigenMode2D(eps, mu, mesh, plane_label, wavelength, num_mode)
            %
            %   Input:
            %     eps - Relative permittivity along each dimension (3D matrix).
            %     mu - Relative permeability along each dimension (3D matrix).
            %     mesh - Grid2D object defining the 2D computational domain.
            %     plane_label - Label indicating the plane ('xy', 'yz', 'zx').
            %     wavelength - Operating wavelength (positive real scalar).
            %     num_mode - Mode numbers to solve for (positive integer array).
            %
            %   Output:
            %     obj - An instance of the `EigenMode2D` class.
            
            arguments
                eps (:,:,3) 
                mu (:,:,3)
                mesh Grid2D
                plane_label {mustBeMember(plane_label, {'xy', 'yz', 'zx'})}
                wavelength (1,1) double {mustBeReal, mustBePositive}
                num_mode (1,:) {mustBePositive, mustBeInteger} = 1
            end

            obj@FDFD2D(mesh);

            % Constant
            eta0 = Constant("eta0").v;
            k0 = 2*pi/wavelength;

            % change all units to m
            obj.mesh = mesh;
            mesh = mesh.changeUnit('m');

            % derivative matrices
            [DE2, DE3] = FDFD.DM2D(mesh, 'E');
            [DH2, DH3] = FDFD.DM2D(mesh, 'M');
            DE2 = DE2/k0;
            DE3 = DE3/k0;
            DH2 = DH2/k0;
            DH3 = DH3/k0;

            % diagonize eps and mu map
            eps1 = sparse(eps(:,:,1)); eps1 = diag(eps1(:));
            eps2 = sparse(eps(:,:,2)); eps2 = diag(eps2(:));
            eps3 = sparse(eps(:,:,3)); eps3 = diag(eps3(:));
            mu1 = sparse(mu(:,:,1)); mu1 = diag(mu1(:));
            mu2 = sparse(mu(:,:,2)); mu2 = diag(mu2(:));
            mu3 = sparse(mu(:,:,3)); mu3 = diag(mu3(:));

            % PQ matrix
            P = [ DE2/eps1*DH3, -(DE2/eps1*DH2+mu3);...
                DE3/eps1*DH3+mu2, -DE3/eps1*DH2 ];
            Q = [ DH2/mu1*DE3, -(DH2/mu1*DE2+eps3);...
                DH3/mu1*DE3+eps2, -DH3/mu1*DE2];

            % eigenmode problem
            eps_core = max(eps,[],"all","ComparisonMethod","real");
            [E,G2] = eigs(P*Q, max(num_mode), -eps_core);
            neff = sqrt(-diag(G2));

            % calculate H field
            H = zeros(size(E));
            for ii = 1:max(num_mode)
                H(:,ii) = 1/(eta0*neff(ii))*Q*E(:,ii);
            end

            % extract and reshape
            e1 = E(1:mesh.axis1.n*mesh.axis2.n,:);
            e2 = E(mesh.axis1.n*mesh.axis2.n+1:end,:);
            h1 = H(1:mesh.axis1.n*mesh.axis2.n,:);
            h2 = H(mesh.axis1.n*mesh.axis2.n+1:end,:);

            E1 = zeros(mesh.axis1.n, mesh.axis2.n, numel(num_mode));
            E2 = zeros(mesh.axis1.n, mesh.axis2.n, numel(num_mode));
            H1 = zeros(mesh.axis1.n, mesh.axis2.n, numel(num_mode));
            H2 = zeros(mesh.axis1.n, mesh.axis2.n, numel(num_mode));
            for ii = 1:numel(num_mode)
                E1(:,:,ii) = reshape(e1(:,num_mode(ii)),[mesh.axis1.n, mesh.axis2.n]);
                E2(:,:,ii) = reshape(e2(:,num_mode(ii)),[mesh.axis1.n, mesh.axis2.n]);
                H1(:,:,ii) = reshape(h1(:,num_mode(ii)),[mesh.axis1.n, mesh.axis2.n]);
                H2(:,:,ii) = reshape(h2(:,num_mode(ii)),[mesh.axis1.n, mesh.axis2.n]);
            end

            % assaign field components
            Ex = zeros(size(E1));
            Ey = zeros(size(E1));
            Ez = zeros(size(E1));
            Hx = zeros(size(E1));
            Hy = zeros(size(E1));
            Hz = zeros(size(E1));
            
            eval(sprintf('E%s = E1;', plane_label(1)));
            eval(sprintf('E%s = E2;', plane_label(2)));
            eval(sprintf('H%s = H1;', plane_label(1)));
            eval(sprintf('H%s = H2;', plane_label(2)));

            obj.lam = wavelength;
            obj.num = num_mode;
            obj.neff = neff(num_mode);
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
                     obj.num;
                     obj.neff};

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength+1, info_name{ii}, value{ii});
            end
        end

        function dispImg(obj, plot_part, varargin)
            arguments
                obj
                plot_part {mustBeMember(plot_part, {'E1', 'E2', 'H1', 'H2'})} = 'E1'
            end
            arguments (Repeating)
                varargin
            end

            field = obj.(plot_part);

            if obj.mesh.axis1.getLength > obj.mesh.axis2.getLength
                row = numel(obj.num);
                col = 1;
            else
                row = ceil(numel(obj.num)/2);
                col = 2;
            end

            % get componet name
            dir = obj.mesh.(['axis',plot_part(2)]).label;
            field_name = [plot_part(1), dir];

            figure;
            sgtitle(field_name);
            for ii = 1:numel(obj.num)
                subplot(row, col, ii)
                plot2D(obj.mesh, field(:,:,ii), varargin{:})
                title(['Mode ', num2str(obj.num(ii)),', neff=',num2str(obj.neff(ii))])
            end
        end

    end

end