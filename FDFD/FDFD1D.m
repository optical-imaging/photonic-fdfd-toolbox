classdef (Abstract) FDFD1D < FDFD
    % FDFD1D Mid-level abstract class of 1D FDFD solvers
    %
    % Description:
    %   The `FDFD1D` class extends the top-level `FDFD` class for use in
    %   one-dimensional FDFD simulations. This class cannot be instantiated directly
    %   and should be subclassed for specific types of 1D simulations.
    %
    % Properties:
    %   mesh - The defined mesh grid along the 1D axis.
    %   lam - The wavelength of the electromagnetic wave in the simulation.
    %   solflag - A flag indicating if the solution has been computed.
    %   Ex, Ey, Ez - Electric field components along x, y, and z dimensions.
    %   Hx, Hy, Hz - Magnetic field components along x, y, and z dimensions.
    %
    % Methods:
    %   FDFD1D(axis) - Constructor to initialize a 1D FDFD solver with an axis.
    %
    % See Also:
    %   FDFD, EigenMode1D

    methods

        function obj = FDFD1D(axis)
            % FDFD1D Constructor
            %   Constructs an instance of the FDFD1D solver with the specified axis.
            %
            %   Syntax:
            %     obj = FDFD1D(axis)
            %
            %   Input:
            %     axis - An axis object of type "Axis" that defines the 1D computational domain.
            %
            %   Output:
            %     obj - An instance of the FDFD1D class.

            arguments
                axis Axis
            end

            obj@FDFD(axis);
            obj.Ex = zeros(axis.n,1);
            obj.Ey = zeros(axis.n,1);
            obj.Ez = zeros(axis.n,1);
            obj.Hx = zeros(axis.n,1);
            obj.Hy = zeros(axis.n,1);
            obj.Hz = zeros(axis.n,1);
        end
    end
    
end