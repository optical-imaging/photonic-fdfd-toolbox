classdef (Abstract) FDFD
    % FDFD Top-level abstract class of FDFD solvers
    %
    % Description:
    %   The `FDFD` class serves as an abstract base class for defining solvers
    %   in the finite-difference frequency-domain (FDFD) method. It cannot be
    %   instantiated directly and should be subclassed for specific types of
    %   FDFD simulations.
    %
    % Properties:
    %   mesh - The defined mesh grid over which the FDFD solution is computed.
    %   lam - The wavelength of the electromagnetic wave in the simulation.
    %   solflag - A flag indicating if the solution has been computed.
    %   Ex, Ey, Ez - Electric field components along x, y, and z dimensions.
    %   Hx, Hy, Hz - Magnetic field components along x, y, and z dimensions.
    %
    % Methods:
    %   FDFD(grid) - Constructor to initialize an FDFD solver with a mesh grid.
    %
    % See Also:
    %   FDFD1D, FDFD2D

    properties (SetAccess = protected)
        mesh
        lam
        solflag
        Ex
        Ey
        Ez
        Hx
        Hy
        Hz
    end

    methods

        function obj = FDFD(grid)
            % FDFD Constructor
            %   Constructs an instance of the FDFD solver with the specified grid.
            %
            %   Syntax:
            %     obj = FDFD(grid)
            %
            %   Input:
            %     grid - A grid object of type "Grid" that defines the computational domain.
            %
            %   Output:
            %     obj - An instance of the FDFD class.
            
            arguments
                grid {mustBeA(grid, "Grid")}
            end

            obj.mesh = grid;
            obj.lam = [];
            obj.solflag = false;
            obj.Ex = [];
            obj.Ey = [];
            obj.Ez = [];
            obj.Hx = [];
            obj.Hy = [];
            obj.Hz = [];
        end
    end

    methods (Static, Access = protected)

        function D = DM1D(axis, field_type)
            switch field_type
                case 'E' % for E components
                    D = spdiags([-ones(axis.n,1) ones(axis.n,1)], [0 1], axis.n, axis.n)/axis.d;
                case 'M' % for H components
                    D = spdiags([-ones(axis.n,1) ones(axis.n,1)], [-1 0], axis.n, axis.n)/axis.d;
            end
        end

        function [D1, D2] = DM2D(mesh, field_type)
            d1 = FDFD.DM1D(mesh.axis1,field_type);
            d2 = FDFD.DM1D(mesh.axis2,field_type);
            D1 = kron(speye(mesh.axis2.n), d1);
            D2 = kron(d2, speye(mesh.axis1.n));
        end

        function [D1, D2, D3] = DM3D(mesh, field_type)
        end
        
    end

end