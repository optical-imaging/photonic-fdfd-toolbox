classdef (Abstract, Hidden) Scattering < FDFD
    %Scattering: Abstract base class for FDFD scattering solvers with external sources
    %
    % Key Properties (inherited from FDFD):
    %   lam         - Wavelength (m)
    %   mesh        - Mesh object (Axis | Grid2D | Grid3D)
    %   bc          - Boundary-condition container
    %   Ex,Ey,Ez    - Solved electric-field components
    %   Hx,Hy,Hz    - Solved magnetic-field components
    %
    % Key Properties (Scattering-specific):
    %   src         - Source object providing incident fields
    %
    % Dependent Properties (inherited):
    %   setflag     - True if lam, mesh, bc, and src are all assigned
    %   solflag     - True if all E/H field components are populated
    %
    % Key Methods (Abstract, Access = protected):
    %   rotateProfile(...)       - Rotate source field profile to match incidence direction
    %   rotatePhase(...)         - Apply phase rotation due to oblique incidence
    %
    % Key Methods:
    %   Scattering(lam)          - Construct scattering solver with wavelength
    %   setSrc(src)              - Attach a configured Source to the solver
    %
    % Key Methods (Access = protected):
    %   checkSet()               - Extend FDFD readiness check to include source

    properties
        src
    end

    methods (Abstract, Access = protected)
        rotateProfile
        rotatePhase
    end

    methods
        % constructor
        function obj = Scattering(wavelength)
            arguments
                wavelength (1,1) double {mustBePositive} % in unit [m]
            end
            obj@FDFD(wavelength);
            obj.src = [];
        end

        % set source
        function setSrc(obj, src)
            arguments
                obj
                src {mustBeA(src, 'Source')}
            end
            if ~src.setflag
                dispError('Scattering:SetSourceFirst');
            end
            if abs(obj.lam-src.lam) > 1e-9
                dispWarning('Scattering:SourceLamNotMatch');
                obj.lam = src.lam;
            end
            obj.src = src;
        end
    end

    methods (Access = protected)
        function temp_status = checkSet(obj)
            temp_status = obj.checkSet@FDFD;
            temp_status = temp_status && ~isempty(obj.src);
        end
    end

end