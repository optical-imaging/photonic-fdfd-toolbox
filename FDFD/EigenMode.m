classdef (Abstract, Hidden) EigenMode < FDFD
    %EigenMode: Abstract base class for FDFD eigenmode solvers
    %
    % Key Properties (inherited from FDFD):
    %   lam         - Wavelength (m)
    %   mesh        - Mesh object (Axis | Grid2D | Grid3D)
    %   bc          - Boundary-condition container
    %   Ex,Ey,Ez    - Solved electric-field eigenmodes
    %   Hx,Hy,Hz    - Solved magnetic-field eigenmodes
    %
    % Key Properties (EigenMode-specific, SetAccess = protected):
    %   modenum     - Mode indices requested/solved (vector of positive integers)
    %   neff        - Effective indices corresponding to solved modes
    %
    % Dependent Properties (inherited):
    %   setflag     - True if lam, mesh, bc, and modenum are assigned
    %   solflag     - True if all E/H field components are populated
    %
    % Key Methods:
    %   EigenMode(lam,modenum)   - Construct eigenmode solver with wavelength and mode indices
    %
    % Key Methods (Access = protected):
    %   checkSet()               - Extend FDFD readiness check to include modenum

    properties (SetAccess = protected)
        modenum
        neff
    end

    methods
        % constructor
        function obj = EigenMode(wavelength, modenum)
            arguments
                wavelength (1,1) double {mustBePositive} % in unit [m]
                modenum (1,:) double {mustBeInteger, mustBePositive}
            end
            obj@FDFD(wavelength);
            obj.modenum = modenum;
            obj.neff = [];
        end
    end

    methods (Access = protected)
        function temp_status = checkSet(obj)
            temp_status = obj.checkSet@FDFD;
            temp_status = temp_status && ~isempty(obj.modenum);
        end
    end

end