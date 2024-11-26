classdef (Abstract) Source < handle
    % Source: FDFD input source for simulations
    %
    % Description:
    %   The `Source` class defines an input source for FDFD (Finite Difference Frequency Domain) simulations.
    %   It provides properties and methods for defining the source characteristics, such as phase, injection direction,
    %   effective refractive index, and electric and magnetic fields.
    %
    % Properties:
    %   setflag - Source set-up flag (boolean to indicate if source is properly set up).
    %   port - The port where the source is injected (Port1D object).
    %   mesh - Defined mesh grid for the source.
    %   lam - Wavelength of the source.
    %   neff - Effective refractive index of the source.
    %   Ex, Ey, Ez - Electric field components along x, y, and z dimensions.
    %   Hx, Hy, Hz - Magnetic field components along x, y, and z dimensions.
    %   phi - Injection phase of the source.
    %   theta - Rotation angle of the source in degrees.
    %
    % Dependent Properties:
    %   dir - Injection direction of the source (derived from the port property).
    %
    % Methods:
    %   Source(port) - Constructor to create an instance of the Source class.
    %   setPhase(phi) - Sets the value of the source phase.
    %   setAngle(theta_in_deg) - Sets the rotation angle of the source in degrees.
    %
    % Abstract Methods:
    %   setEigenMode - Abstract method to set an eigenmode as the source.
    %   setPlaneWave - Abstract method to set a plane wave as the source.
    %   setGaussianBeam - Abstract method to set a Gaussian beam as the source.
    %
    % See also:
    %   Port1D, Device, Device2D, Grid, Axis, Scattering2D

    properties
        setflag
        port
        mesh
        lam
        neff
        Ex
        Ey
        Ez
        Hx
        Hy
        Hz
        phi
        theta
    end

    properties (Dependent)
        dir
    end

    methods (Abstract)
        setEigenMode
        setPlaneWave
        setGaussianBeam
    end

    methods

        function obj = Source(port)
            % Source constructor
            %   Constructs an instance of the Source class with the specified port.
            %
            %   Syntax:
            %     obj = Source(port)
            %
            %   Input:
            %     port - The port where the source is injected (Port1D object).
            %
            %   Output:
            %     obj - An instance of the Source class.

            obj.setflag = false;
            obj.port = port;
            obj.mesh = [];
            obj.lam = [];
            obj.neff = [];
            obj.Ex = [];
            obj.Ey = [];
            obj.Ez = [];
            obj.Hx = [];
            obj.Hy = [];
            obj.Hz = [];
            obj.phi = 0;
            obj.theta = 0;
        end

        function value = get.dir(obj)
            value = obj.port.dir;
        end

        function obj = setPhase(obj, phi)
            % setPhase - Sets the injection phase of the source
            %
            %   Syntax:
            %     obj = setPhase(phi)
            %
            %   Input:
            %     phi - Phase value to set (real number).

            arguments
                obj 
                phi (1,1) {mustBeReal}
            end

            obj.phi = phi;
        end

        function obj = setAngle(obj, theta_in_deg)
            % setAngle - Sets the rotation angle of the source in degrees
            %
            %   Syntax:
            %     obj = setAngle(theta_in_deg)
            %
            %   Input:
            %     theta_in_deg - Rotation angle in degrees (must be in range -180 to 180).

            arguments
                obj 
                theta_in_deg (1,1) {mustBeInRange(theta_in_deg, -180, 180)}
            end

            obj.theta = theta_in_deg;
        end

    end
end