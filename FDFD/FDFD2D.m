classdef (Abstract) FDFD2D < FDFD
     % FDFD2D: Mid-level abstract class of 2D FDFD solvers
    %
    % Description:
    %   The `FDFD2D` class extends the top-level `FDFD` class for use in
    %   two-dimensional FDFD simulations. This class cannot be instantiated directly
    %   and should be subclassed for specific types of 2D simulations.
    %
    % Properties:
    %   mesh - The defined mesh grid for the 2D domain.
    %   lam - The wavelength of the electromagnetic wave in the simulation.
    %   solflag - A flag indicating if the solution has been computed.
    %   Ex, Ey, Ez - Electric field components along x, y, and z dimensions.
    %   Hx, Hy, Hz - Magnetic field components along x, y, and z dimensions.
    %
    % Methods:
    %   FDFD2D(grid) - Constructor to initialize a 2D FDFD solver with a grid.
    %   rotateProfile(source_field, axis_source, axis_prop, theta) - Static method to rotate a source field profile.
    %   rotatePhase(beta, AXIS_ds_source, AXIS_ds_prop, port_position, theta, phi) - Static method to rotate a phase profile.
    %
    % Example:
    %   % Define a grid for the 2D simulation
    %   grid = Grid2D(Axis('x', 0, 1, 100), Axis('y', 0, 1, 100));
    %
    %   % Create a subclass instance of FDFD2D (e.g., a specific 2D solver)
    %   solver = MyFDFD2DSolver(grid);
    %
    % Notes:
    %   - This class is abstract and should be subclassed for specific 2D simulation purposes.
    %   - Ensure that the `Grid2D` object is properly defined before using it.
    %
    % See Also:
    %   FDFD, Grid, Grid2D, EigenMode2D, Scattering2D

    methods

        function obj = FDFD2D(grid)
            % FDFD2D: Constructor to initialize a 2D FDFD solver with a grid
            %
            %   Syntax:
            %     obj = FDFD2D(grid)
            %
            %   Input:
            %     grid - A `Grid2D` object that defines the 2D computational domain.
            %
            %   Output:
            %     obj - An instance of the `FDFD2D` class.
            
            arguments
                grid Grid2D
            end

            obj@FDFD(grid);
            obj.Ex = zeros(grid.axis1.n, grid.axis2.n);
            obj.Ey = zeros(grid.axis1.n, grid.axis2.n);
            obj.Ez = zeros(grid.axis1.n, grid.axis2.n);
            obj.Hx = zeros(grid.axis1.n, grid.axis2.n);
            obj.Hy = zeros(grid.axis1.n, grid.axis2.n);
            obj.Hz = zeros(grid.axis1.n, grid.axis2.n);
        end

    end

    methods (Static, Access = protected)

        function profile_rot = rotateProfile(source_field, axis_source, axis_prop, theta)
            arguments
                source_field (1,:)
                axis_source
                axis_prop
                theta
            end

            % rotate grid
            theta = theta/180*pi;
            [PROP, SOURCE] = ndgrid(axis_prop.v, axis_source.v);
            [TH, R] = cart2pol(PROP, SOURCE);
            [~, SOURCE] = pol2cart(TH+theta, R);

            % interpolate source
            profile_interp = griddedInterpolant(axis_source.v, source_field,...
                'linear', 'none');
            profile_rot = profile_interp(SOURCE);
            profile_rot(isnan(profile_rot)) = 0;
        end

        function phase2_rot = rotatePhase(beta, AXIS_ds_source, AXIS_ds_prop, port_position, theta, phi)
            % rotate grid
            theta = theta/180*pi;
            PROP = AXIS_ds_prop-port_position;
            SOURCE = AXIS_ds_source;
            [TH, R] = cart2pol(PROP, SOURCE);
            [PROP, ~] = pol2cart(TH+theta, R);

            phase2_rot = exp(-1i*beta*PROP+1i*phi);
        end

    end

end