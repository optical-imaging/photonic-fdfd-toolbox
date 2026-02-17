classdef (Abstract, Hidden) Source < handle
    %Source: Abstract base class for defining excitation fields for FDFD simulations
    %
    % Key Properties:
    %   lam         - Source wavelength (m)
    %   port        - Port object defining injection/projection location (provides dir)
    %   mesh        - Source mesh extracted from the simulation mesh around the port
    %   neff        - Effective index of the source mode (sign corrected by dir)
    %   phi         - Global phase offset (rad)
    %   theta       - Incident angle (deg), used by plane-wave / Gaussian-beam sources
    %
    %   Ex,Ey,Ez    - Electric-field components of the source (port-mesh sampled)
    %   Hx,Hy,Hz    - Magnetic-field components of the source (port-mesh sampled)
    %
    % Dependent Properties:
    %   dir         - Injection direction, forwarded from port.dir
    %   setflag     - True if mesh/neff and all E/H components are populated
    %
    % Key Methods (Abstract):
    %   setMesh(mesh)                 - Bind a simulation mesh and build the port mesh
    %   setEigenMode(device, ...)     - Configure source using an eigenmode on the port cut
    %   setPlaneWave(BG_material, ...) - Configure source as a plane wave in background
    %   setGaussianBeam(BG_material, ...) - Configure source as a Gaussian beam in background
    %
    % Key Methods (Abstract, Static, Access = protected):
    %   getPortMesh(port, mesh)       - Extract a 1D/2D port-aligned mesh from the full mesh
    %
    % Key Methods (Abstract, Access = {?Scattering, ?Source}):
    %   extendField(axis)             - Embed/extend port fields onto a larger axis/grid
    %
    % Key Methods:
    %   Source(lam, port)             - Construct source container and initialize defaults
    %   get.dir()                     - Return port.dir
    %   get.setflag()                 - Check whether fields and parameters are ready
    %   setPhase(phi)                 - Set phase offset
    %   setAngle(theta_deg)           - Set incident angle in degrees
    %   setAmplitude(amp, seq, 'E'|'H') - Scale selected source instance(s) by field magnitude
    %
    % Key Methods (Access = protected):
    %   assignEigenMode(eigenmode1d)  - Copy fields from EigenMode1D and correct neff sign
    %   correctNeff(neff)             - Flip neff sign based on injection direction
    %
    % Key Methods (Static, Access = private):
    %   indexHelper(obj, ii)          - Return indexing cell for Source1D/Source2D field storage

    properties
        lam
        port
        mesh
        neff
        phi
        theta
        Ex
        Ey
        Ez
        Hx
        Hy
        Hz
    end

    properties (Dependent)
        dir
        setflag
    end

    methods (Abstract)
        setMesh
        setEigenMode
        setPlaneWave
        setGaussianBeam
    end

    methods (Abstract, Static, Access = protected)
        getPortMesh
    end

    methods (Abstract, Access = {?Scattering, ?Source})
        extendField
    end

    methods
        % constructor
        function obj = Source(wavelength, port)
            arguments
                wavelength (1,1) double {mustBePositive}
                port
            end

            obj.lam = wavelength;
            obj.port = port;
            obj.mesh = [];
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

        % dependent
        function val = get.dir(obj)
            val = obj.port.dir;
        end

        function val = get.setflag(obj)
            val = ~isempty(obj.mesh)  && ...
                ~isempty(obj.neff)  && ...
                ~isempty(obj.Ex)    && ...
                ~isempty(obj.Ey)    && ...
                ~isempty(obj.Ez)    && ...
                ~isempty(obj.Hx)    && ...
                ~isempty(obj.Hy)    && ...
                ~isempty(obj.Hz);
        end

        % manipulate
        function setPhase(obj, phi)
            arguments
                obj
                phi (1,1) {mustBeReal}
            end
            obj.phi = phi;
        end

        function setAngle(obj, theta_in_deg)
            arguments
                obj
                theta_in_deg (1,1) {mustBeInRange(theta_in_deg, -180, 180)}
            end
            obj.theta = theta_in_deg;
        end

        function setAmplitude(obj, amp, src_seq, fields)
            arguments
                obj
                amp (1,1) double {mustBeGreaterThan(amp, 0)}
                src_seq (1,:) double {mustBeInteger, mustBePositive} = 1
                fields char {mustBeMember(fields, {'E', 'H'})} = 'E'
            end

            for ii = src_seq
                idx = obj.indexHelper(ii);  % cell, like {':', ii} or {':',':',ii}

                switch fields
                    case 'E'
                        mag = sqrt( abs(obj.Ex(idx{:})).^2 + ...
                            abs(obj.Ey(idx{:})).^2 + ...
                            abs(obj.Ez(idx{:})).^2 );
                    case 'H'
                        mag = sqrt( abs(obj.Hx(idx{:})).^2 + ...
                            abs(obj.Hy(idx{:})).^2 + ...
                            abs(obj.Hz(idx{:})).^2 );
                end

                mag0 = max(mag(:));
                if mag0 == 0
                    error('Source:setAmplitude:ZeroMagnitude', ...
                        'Cannot scale source #%d because reference magnitude is zero.', ii);
                end

                scale_factor = amp / mag0;

                % scale all components for that source instance
                fieldNames = {'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'};
                for f = 1:numel(fieldNames)
                    fname = fieldNames{f};
                    obj.(fname)(idx{:}) = obj.(fname)(idx{:}) * scale_factor;
                end
            end
        end
    end

    methods (Access = protected)
        function assignEigenMode(obj, eigenmode)
            arguments
                obj
                eigenmode EigenMode1D
            end

            obj.correctNeff(eigenmode.neff);
            obj.Ex   = eigenmode.Ex;
            obj.Ey   = eigenmode.Ey;
            obj.Ez   = eigenmode.Ez;
            obj.Hx   = eigenmode.Hx;
            obj.Hy   = eigenmode.Hy;
            obj.Hz   = eigenmode.Hz;
        end

        function correctNeff(obj, neff)
            switch obj.dir(1)
                case '+'
                    obj.neff = neff;
                case '-'
                    obj.neff = -neff;
            end
        end
    end

    methods (Access = private)
        function idx = indexHelper(obj, ii)
            switch class(obj)
                case 'Source1D'
                    idx = {':', ii};
                case 'Source2D'
                    idx = {':', ':', ii};
            end
        end
    end
end
