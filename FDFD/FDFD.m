classdef (Abstract, Hidden) FDFD < handle
    %FDFD: Abstract base class for frequency-domain solvers operating on a mesh with BCs
    %
    % Key Properties (SetAccess = protected):
    %   lam         - Wavelength (m)
    %   mesh        - Mesh object (Axis | Grid2D | Grid3D)
    %   bc          - Boundary-condition container (BC1D | BC2D | BC3D)
    %   Ex,Ey,Ez    - Electric-field components of the solved fields
    %   Hx,Hy,Hz    - Magnetic-field components of the solved fields
    %
    % Dependent Properties:
    %   setflag     - True if lam, mesh, and bc are all assigned (solver ready)
    %   solflag     - True if all E/H field components are populated (solution available)
    %
    % Key Methods (Abstract):
    %   solveFDFD(...)            - Solve frequency-domain Maxwell equations
    %
    % Key Methods:
    %   FDFD(lam)                 - Construct solver container with wavelength
    %   get.setflag()             - Return solver readiness flag
    %   get.solflag()             - Return solver solution flag
    %   setMesh(mesh)             - Assign mesh to solver
    %   setBC(bc)                 - Assign boundary conditions to solver
    %
    % Key Methods (Access = protected):
    %   checkSet()                - Internal readiness check used by setflag
    %
    % Key Methods (Static, Access = protected):
    %   DM1D(axis,type)           - 1D finite-difference derivative (Dirichlet BC)
    %   DM2D(grid2d,type)         - 2D finite-difference derivatives (Dirichlet BC)
    %   DM3D(grid3d)              - 3D finite-difference derivatives (Dirichlet BC)
    %   DM1D_PBC(axis,type,kinc)  - 1D finite-difference derivative (Periodic BC)
    %   DM2D_PBC(grid2d,type,kinc)- 2D finite-difference derivatives (Periodic BC)
    %   DM3D_PBC(grid3d)          - 3D finite-difference derivatives (Periodic BC)

    properties (SetAccess = protected)
        lam
        mesh
        bc
        Ex
        Ey
        Ez
        Hx
        Hy
        Hz
    end

    properties (Dependent)
        setflag
        solflag
    end

    methods (Abstract)
        setMesh
        setBC
        solveFDFD
    end

    methods
        % constructor
        function obj = FDFD(wavelength)
            arguments
                wavelength (1,1) double {mustBePositive} % in unit [m]
            end

            obj.lam = wavelength;
            obj.mesh = [];
            obj.bc = [];
            obj.Ex = []; obj.Ey = []; obj.Ez = [];
            obj.Hx = []; obj.Hy = []; obj.Hz = [];
        end

        % dependent
        function val =  get.setflag(obj)
            val = obj.checkSet;
        end

        function val = get.solflag(obj)
            val = ~isempty(obj.Ex) ...
                  && ~isempty(obj.Ey) ...
                  && ~isempty(obj.Ez) ...
                  && ~isempty(obj.Hx) ...
                  && ~isempty(obj.Hy) ...
                  && ~isempty(obj.Hz);
        end
    end

    methods (Access = protected)
        function temp_status = checkSet(obj)
            temp_status = ~isempty(obj.lam)...
                          && ~isempty(obj.mesh)...
                          && ~isempty(obj.bc);
        end
    end

    methods (Static, Access = protected)
        % 1D finite difference matrix (Dirichlet)
        function D = DM1D(axis, field_type)
            switch field_type
                case 'E'
                    D = spdiags([-ones(axis.n,1), ones(axis.n,1)], [0 1], axis.n, axis.n) / axis.d;
                case 'M'
                    D = spdiags([-ones(axis.n,1), ones(axis.n,1)], [-1 0], axis.n, axis.n) / axis.d;
            end
        end

        % 2D finite difference matrices (Dirichlet)
        function [D1, D2] = DM2D(mesh, field_type)
            d1 = FDFD.DM1D(mesh.axis1, field_type);
            d2 = FDFD.DM1D(mesh.axis2, field_type);
            D1 = kron(speye(mesh.axis2.n), d1);
            D2 = kron(d2, speye(mesh.axis1.n));
        end

        % 3D finite difference matrices (Dirichlet)
        function [Dx, Dy, Dz] = DM3D(mesh)
        end

        % 1D finite difference matrix (Periodic BC)
        function D = DM1D_PBC(axis, field_type, k_inc)
            L = axis.n * axis.d;
            switch field_type
                case 'E'
                    D = spdiags([-ones(axis.n,1), ones(axis.n,1)], [0 1], axis.n, axis.n) / axis.d;
                    D(end,1) = exp(-1i * k_inc * L) / axis.d;
                case 'M'
                    D = spdiags([-ones(axis.n,1), ones(axis.n,1)], [-1 0], axis.n, axis.n) / axis.d;
                    D(1,end) = -exp(1i * k_inc * L) / axis.d;
            end
        end

        % 2D finite difference matrices (Periodic BC)
        function [D1, D2] = DM2D_PBC(mesh, field_type, k_inc)
            d1 = FDFD.DM1D_PBC(mesh.axis1, field_type, k_inc(1));
            d2 = FDFD.DM1D_PBC(mesh.axis2, field_type, k_inc(2));
            D1 = kron(speye(mesh.axis2.n), d1);
            D2 = kron(d2, speye(mesh.axis1.n));
        end

        % 3D finite difference matrices (Periodic BC)
        function [Dx, Dy, Dz] = DM3D_PBC(mesh)
        end
    end
end
