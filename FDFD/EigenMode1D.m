classdef (Hidden) EigenMode1D < FDFD1D
    % EigenMode1D: 1D eigenmode solver for electromagnetic analysis
    %
    % Description:
    %   The `EigenMode1D` class is used for computing the eigenmodes of one-dimensional electromagnetic structures.
    %   It provides capabilities to solve for the effective refractive index and the field distributions (electric and magnetic fields).
    %   This class inherits from `FDFD1D` and utilizes finite-difference frequency-domain (FDFD) techniques for calculations.
    %
    % Properties:
    %   mesh - The defined mesh grid along the 1D axis.
    %   lam - The wavelength of the electromagnetic wave in the simulation (in meters).
    %   solflag - A flag indicating if the solution has been computed.
    %   Ex, Ey, Ez - Electric field components along x, y, and z dimensions.
    %   Hx, Hy, Hz - Magnetic field components along x, y, and z dimensions.
    %   pol - Polarization type ('TE' or 'TM').
    %   num - Mode number to solve for.
    %   neff - Effective refractive index of the computed mode.
    %
    % Methods:
    %   EigenMode1D(eps, mu, axis, axis_label, wavelength, Name, Value) - Constructor to create a 1D eigenmode solver object.
    %   printInfo() - Displays information about the computed mode.
    %
    % Example:
    %   % Define material properties, mesh grid, and wavelength
    %   eps = [2.5, 2.5, 2.5];
    %   mu = [1.0, 1.0, 1.0];
    %   axis = Axis('x', 0, 1, 100);
    %   axis_label = 'x';
    %   wavelength = 1.55e-6;
    %
    %   % Create an EigenMode1D solver object and compute the TE mode
    %   emode = EigenMode1D(eps, mu, axis, axis_label, wavelength, 'ModeType', 'TE', 'ModeNum', 1);
    %   emode.printInfo();
    %
    % Notes:
    %   - This class is marked as "Hidden", which means it is not intended to be directly accessed by users.
    %   - Ensure that the axis object is defined properly, with the units set to meters.
    %
    % See Also:
    %   FDFD, FDFD1D, Grid, Axis

    properties (SetAccess = private)
        pol
        num
        neff
    end

    methods

        function obj = EigenMode1D(eps, mu, axis, axis_label, wavelength, varargin)
            % EigenMode1D: Constructor to create a 1D eigenmode solver object
            %
            %   Syntax:
            %     obj = EigenMode1D(eps, mu, axis, axis_label, wavelength, Name, Value)
            %
            %   Input:
            %     eps - Relative permittivity along each dimension (matrix).
            %     mu - Relative permeability along each dimension (matrix).
            %     axis - Axis object containing the mesh grid.
            %     axis_label - Label indicating the direction of propagation ('x', 'y', or 'z').
            %     wavelength - Operating wavelength (positive real scalar).
            %
            %   Name-Value Arguments:
            %     'ModeType' - Type of mode ('TE' or 'TM'). Default is 'TE'.
            %     'ModeNum' - Mode number to solve for. Default is 1.
            
            arguments
                eps (:,3) 
                mu (:,3)
                axis Axis
                axis_label {mustBeMember(axis_label, {'x', 'y', 'z'})}
                wavelength (1,1) double {mustBeReal, mustBePositive}
            end
            arguments (Repeating)
                varargin
            end

            obj@FDFD1D(axis);

            % Set up inputParser to handle name-value pairs
            p = inputParser;
            addParameter(p, 'ModeType', 'TE', @(x) ismember(x, {'TE', 'TM'})); % Default polarization mode is TE
            addParameter(p, 'ModeNum', 1, @(x) isPositiveInt(x));
            parse(p, varargin{:});
            pol_mode = p.Results.ModeType;
            n_mode = sort(p.Results.ModeNum);

            % Constant
            eta0 = Constant("eta0").v;
            k0 = 2*pi/wavelength;

            % change axis unit to m
            obj.mesh = axis; % remain original unit before changing
            axis = axis.changeUnit('m');

            % derivative matrices
            DE2 = FDFD.DM1D(axis, 'E')/k0;
            DH2 = FDFD.DM1D(axis, 'M')/k0;

            % solve eigenmode
            switch pol_mode
                case 'TE'
                    A = -(DH2*diag(sparse(mu(:,1)))*DE2+diag(sparse(eps(:,3))));
                    [E3, G2] = eigs(A, diag(sparse(1./mu(:,2))), max(n_mode), -max(eps,[],"all"));
                    neff = sqrt(-diag(G2));

                    H1 = diag(sparse(1./mu(:,1)))*DE2*E3;
                    H2 = sqrt(diag(G2)).'.*(diag(sparse(1./mu(:,2)))*E3);

                    [~, ind] = sort(real(neff),'descend');
                    e3 = E3(:,ind);
                    h1 = H1(:,ind)./(-1i*eta0);
                    h2 = H2(:,ind)./(-1i*eta0); % out of phase for Ez and Hy, because wave is propagating towards +x
                    neff = neff(ind);

                    e3 = e3(:,n_mode);
                    h1 = h1(:,n_mode);
                    h2 = h2(:,n_mode);

                    e1 = zeros(size(eps,1), numel(n_mode));
                    e2 = zeros(size(eps,1), numel(n_mode));
                    h3 = zeros(size(mu,1), numel(n_mode));
                case 'TM'
                    A = -(DE2*diag(sparse(eps(:,1)))*DH2+diag(sparse(mu(:,3))));
                    [H3, G2] = eigs(A, diag(sparse(1./eps(:,2))), max(n_mode), -max(eps,[],"all"));
                    neff = sqrt(-diag(G2));

                    E1 = diag(sparse(1./eps(:,1)))*DH2*H3;
                    E2 = sqrt(diag(G2)).'.*(diag(sparse(1./eps(:,2)))*H3);

                    [~, ind] = sort(real(neff),'descend');
                    h3 = H3(:,ind)./(-1i*eta0)/1i; % fix phase cause in EY
                    e1 = E1(:,ind)/1i;
                    e2 = E2(:,ind)/1i;
                    neff = neff(ind);

                    h3 = h3(:,n_mode);
                    e1 = e1(:,n_mode);
                    e2 = e2(:,n_mode);

                    e3 = zeros(size(eps,1), numel(n_mode));
                    h1 = zeros(size(mu,1), numel(n_mode));
                    h2 = zeros(size(mu,1), numel(n_mode));
            end

            [E, H] = EigenMode1D.assignField(axis_label, e1, e2, e3, h1, h2, h3);

            obj.lam = wavelength;
            obj.solflag = true;
            obj.pol = pol_mode;
            obj.num = n_mode;
            obj.neff = neff(n_mode);
            obj.Ex = E.Ex;
            obj.Ey = E.Ey;
            obj.Ez = E.Ez;
            obj.Hx = H.Hx;
            obj.Hy = H.Hy;
            obj.Hz = H.Hz;
        end

        function printInfo(obj)
            info_name = {'Object name';
                         'Wavelength';
                         'Polarization mode';
                         'Mode number';
                         'neff'};
            value = {inputname(1);
                     [num2str(obj.lam), ' m'];
                     obj.pol;
                     obj.num;
                     obj.neff};

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength+1, info_name{ii}, value{ii});
            end
        end

    end

    methods (Static, Access = private)

        function [E, H] = assignField(axis_label, e1, e2, e3, h1, h2, h3)
            perm = struct('x', [2,3,1], 'y', [1,2,3], 'z', [3,1,2]);
            p = perm.(axis_label);

            E.Ex = eval(sprintf('e%d', p(1))); H.Hx = eval(sprintf('h%d', p(1)));
            E.Ey = eval(sprintf('e%d', p(2))); H.Hy = eval(sprintf('h%d', p(2)));
            E.Ez = eval(sprintf('e%d', p(3))); H.Hz = eval(sprintf('h%d', p(3)));
        end

    end



end