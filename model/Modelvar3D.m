classdef Modelvar3D < Model
    %Modelvar3D: Variational 3D workflow model that stacks LayerDevice (z-separated) on Grid3D and solves with varFDFD
    %
    % Key Properties (inherited from Model, SetAccess = protected):
    %   stage       - Workflow stage/state
    %   lam         - Wavelength (m)
    %   bgmat       - Background material (Material)
    %   devicelist  - LayerDevice component list (indexed by seq)
    %   device      - Assembled LayerDevice used for simulation (layer stack)
    %   mesh        - Grid3D mesh (x,y,z)
    %   solver      - Solver selection (fixed as 'var' here)
    %   cfg         - Configuration struct (bc, pol, port, source, etc.)
    %   solution    - varFDFD solver object (stores epseff, ref, and solved 3D fields)
    %
    % Key Methods (Initialization / reset):
    %   Modelvar3D()                      - Constructor (called only by Model.initializeSingleObj)
    %   resetProperties()                  - Reset Semi-3D state and allocate empty LayerDevice containers
    %   initialize()                       - Return singleton Modelvar3D instance
    %
    % Key Methods (Device assembly):
    %   assembleDevice()                   - Assemble LayerDevice list (calls Model.assembleDevice) and
    %                                       check z-overlap via device.layerOverlap()
    %
    % Key Methods (Meshing / coordinates):
    %   setMesh(num_point,step_size,...)   - Create Grid3D mesh (signature enforces 1×3 inputs)
    %   getCoordinate()                    - Return x,y,z coordinate vectors from Grid3D mesh
    %   addDeviceBit(eps_bit,mu_bit,step)  - Define LayerDevice directly from 3D bitmap arrays
    %
    % Key Methods (Solver setup):
    %   setSolver(Name,Value,...)          - Configure varFDFD (the only solver):
    %                                       - Default BC: PML on x and y (BC2D on 'xy')
    %                                       - Parameters: 'Pol' ('TE'|'TM'), 'RefPoint' ([x y])
    %                                       - Calls: solution.setMesh(Grid3D), solution.setBC(BC2D),
    %                                         solution.setReference(LayerDevice,RefPoint)
    %
    % Key Methods (Boundary condition):
    %   setBoundary(axis_label,type,opt)   - Update BC2D on x or y; for PML, opt sets thickness
    %
    % Key Methods (Ports and sources for varFDFD driving field):
    %   addPort(direction,position,span)   - Create Port1D on 'xy' plane and clear existing source
    %   addSource(source_type,opt)         - Build Source1D on the effective 2D device (epseff):
    %                                       - Creates Device2D "Alpha" from solution.epseff on mesh.getPlane('xy')
    %                                       - Source types: 'EigenMode' (opt modes), 'GaussianBeam' (opt waist),
    %                                         'PlaneWave'
    %                                       - Attaches to varFDFD via solution.setSrc(Source1D)
    %
    % Key Methods (Run and output):
    %   solveModel()                       - Run varFDFD solve (solution.solveFDFD) after port+source setup
    %   exportResults()                    - Return solver object (varFDFD) as results container
    %   dispInfo()                         - Placeholder (not implemented)
    %
    % Private helper (currently unused in main flow):
    %   initvarFDFD(Name,Value,...)        - Older initialization path for varFDFD (kept for internal refactoring)

    methods (Access = private)
        % constructor
        function obj = Modelvar3D()
            obj.resetProperties();
        end

        function resetProperties(obj)
            obj.resetPartialProperties;
            obj.devicelist = LayerDevice.empty(0,1);
            obj.device = LayerDevice.empty(0,1);
        end
    end

    methods (Static)
        function singleObj = initialize()
            persistent uniqueInstance

            if isempty(uniqueInstance) || ~isvalid(uniqueInstance)
                disp("Modelvar3D: Object Initialized");
            else
                disp("Modelvar3D: Object Reset and Initialized");
            end

            uniqueInstance = Modelvar3D();
            singleObj = uniqueInstance;
        end
    end

    methods
        % assemble device and check overlapping
        function assembleDevice(obj)
            obj.assembleDevice@Model;

            % check overlapping
            [isoverlap, overlap_seq] = obj.device.layerOverlap;
            if isoverlap
                dispMessage('LayerDevice:LayerOverlap', overlap_seq);
            end
        end

        % meshing
        function setMesh(obj, num_point, step_size)
            arguments
                obj
                num_point (1,3) double {mustBePositive, mustBeInteger}
                step_size (1,3) double {mustBeReal, mustBePositive}
            end
            obj.checkStage([2,3,4], 'Model:SetMesh');

            obj.mesh = Grid3D(num_point, step_size);
            obj.stage = 3;
        end

        function [x_val, y_val, z_val] = getCoordinate(obj)
            if isempty(obj.mesh)
                dispError('Model:NoMesh');
            end
            x_val = obj.mesh.axisx.v;
            y_val = obj.mesh.axisy.v;
            z_val = obj.mesh.axisz.v;
        end

        % Directly set bitmap
        function addDeviceBit(obj, eps_bit, mu_bit, step_size)
            arguments
                obj
                eps_bit (:,:,:) double {mustBeRealAndGE1}
                mu_bit (:,:,:) double {mustBeRealAndGE1}
                step_size (1,3) double {mustBeReal, mustBePositive}
            end
            
            obj.checkStage(1, 'Model:AddDevice');

            obj.device = LayerDevice(1);
            obj.device.setBitmap(eps_bit, mu_bit, step_size);

            obj.mesh = obj.device.mesh;

            obj.stage = 4;
        end

        % Solver setup
        function setSolver(obj, varargin)
            arguments
                obj
            end
            arguments(Repeating)
                varargin
            end

            obj.checkStage([4,5,6,7,8], 'Modelvar3D:SetSolver');

            if isempty(obj.lam)
                dispError('Model:SetWavelength');
            end

            % varFDFD is the only solver
            obj.solver = 'var';
            obj.cfg = struct();   % reset cfg

            % PML is default boundary condition for Scattering solver
            bc1 = BC1D('x', 'PML'); % default PML has 8 layers
            bc2 = BC1D('y', 'PML');
            obj.cfg.bc = BC2D(bc1, bc2);

            % initialize and configure varFDFD
            p = inputParser;
            addParameter(p, 'Pol', 'TE', @(x) ismember(x, {'TE','TM'}));
            addParameter(p, 'RefPoint', [0 0], @(x) isnumeric(x) && numel(x)==2);
            parse(p, varargin{:});

            pol  = p.Results.Pol;
            rpin = p.Results.RefPoint;

            obj.solution = varFDFD(obj.lam, pol);
            obj.solution.setMesh(obj.mesh);
            obj.solution.setBC(obj.cfg.bc);
            obj.solution.setReference(obj.device, rpin);

            obj.cfg.pol = obj.solution.pol;

            obj.stage = 5;
        end

        % Boundary condition
        function setBoundary(obj, axis_label, BC_type, opt)
            arguments
                obj
                axis_label {mustBeMember(axis_label, {'x','y'})}
                BC_type {mustBeMember(BC_type, {'Perfect', 'PML', 'Periodic'})}
                opt = [] % only for setting number of layers in PML
            end

            obj.checkStage([5,6,7,8], 'Modelvar3D:SetBC');

            if strcmp(BC_type, "Perfect")
                dispWarning('Model:PerfectForScattering'); % Perfect BC can cause error in Scattering problem.
            end

            % set new BC
            obj.cfg.bc.changeType(axis_label, BC_type);
            if ~isempty(opt) && strcmp(BC_type, "PML")
                obj.cfg.bc.setPMLThickness(axis_label, opt);
            end

            obj.solution.setBC(obj.cfg.bc); % update BC in Scattering2D if changed

            if obj.stage == 8 % model has been solved
                obj.stage = 7; % back to before solving
            end
        end

        % Source/Port
        function addPort(obj, direction, position, span_range)
            arguments
                obj
                direction char {mustBeMember(direction, {'+x', '-x', '+y', '-y'})}
                position double {mustBeReal}
                span_range (1,2) double {mustBeReal}
            end

            obj.checkStage([5,6,7,8], 'Model:AddPort');

            % Port has be in x-y plane
            if ~contains('xy', direction(2))
                dispError('Port1D:DirectionWrong', 'x', 'y');
            end

            obj.cfg.port = Port1D('xy', direction, position, span_range);
            obj.cfg.source = []; % clear source if port has been changed

            obj.stage = 6;
        end

        function addSource(obj, num_mode)
            arguments
                obj
                num_mode (1,:) double {mustBeInteger, mustBePositive} = 1
            end

            obj.checkStage([6,7,8], 'Modelvar3D:AddSource');

            % create fake 2d device alpha_eff
            gridxy = obj.mesh.getPlane('xy');
            meshxy_step = [obj.mesh.axisx.d, obj.mesh.axisy.d];

            Alpha = Device2D(1);
            eps_eff = obj.solution.epseff;
            mu_eff = ones(size(eps_eff));
            Alpha.setBitmap(eps_eff,mu_eff,meshxy_step,gridxy.label);
            Alpha.shiftMesh('x', 'min', gridxy.axis1.v(1)); % align mesh
            Alpha.shiftMesh('y', 'min', gridxy.axis2.v(1));

            % eigenmode source based on effective index
            s = Source1D(obj.lam, obj.cfg.port, obj.cfg.pol);
            s.setMesh(gridxy);
            s.setEigenMode(Alpha, num_mode); 

            obj.cfg.source = s;
            obj.solution.setSrc(s);

            obj.stage = 7;
        end

        % solve
        function solveModel(obj)
            obj.checkStage(7, 'Modelvar3D:SetBeforeSolving');
            obj.solution.solveFDFD;
            obj.stage = 8;
        end

        % export
        % dispInfo  Show a summary of the current Modelvar3D configuration (developer debugging)
        function dispInfo(obj)
            fprintf('=== Modelvar3D Information ===\n');

            % --- General ---
            if isempty(obj.lam)
                fprintf('Wavelength   : Not set\n');
            else
                fprintf('Wavelength   : %.3e m\n', obj.lam);
            end

            if isempty(obj.bgmat)
                fprintf('Background   : Not set\n');
            else
                fprintf('Background   : eps_r=%.3f, mu_r=%.3f\n', obj.bgmat.eps, obj.bgmat.mu);
            end

            fprintf('Solver       : %s\n', string(obj.solver));

            % --- Mesh ---
            if isempty(obj.mesh)
                fprintf('Mesh         : Not set\n');
            else
                fprintf('Mesh (x)     : N=%d, d=%.3e m\n', obj.mesh.axisx.n, obj.mesh.axisx.d);
                fprintf('Mesh (y)     : N=%d, d=%.3e m\n', obj.mesh.axisy.n, obj.mesh.axisy.d);
                fprintf('Mesh (z)     : N=%d, d=%.3e m\n', obj.mesh.axisz.n, obj.mesh.axisz.d);
            end

            % --- Boundary conditions (xy plane) ---
            if isempty(obj.cfg) || ~isfield(obj.cfg,'bc') || isempty(obj.cfg.bc)
                fprintf('Boundary Cond: Not set\n');
            else
                bcx = obj.cfg.bc.bc1;   % fixed: x
                bcy = obj.cfg.bc.bc2;   % fixed: y

                if strcmp(bcx.type,'PML')
                    fprintf('BC (x)       : %s, thickness=%d\n', bcx.type, bcx.pmlnum);
                else
                    fprintf('BC (x)       : %s\n', bcx.type);
                end

                if strcmp(bcy.type,'PML')
                    fprintf('BC (y)       : %s, thickness=%d\n', bcy.type, bcy.pmlnum);
                else
                    fprintf('BC (y)       : %s\n', bcy.type);
                end
            end


            % --- varFDFD reference ---
            if ~isempty(obj.solution) && isprop(obj.solution,'ref') && ~isempty(obj.solution.ref)
                r = obj.solution.ref;
                if isfield(r,'pin')
                    fprintf('Ref point    : [%.3e, %.3e] m\n', r.pin(1), r.pin(2));
                else
                    fprintf('Ref point    : (not set)\n');
                end
            else
                fprintf('Reference    : Not set\n');
            end

            % --- Polarization ---
            if isfield(obj.cfg,'pol') && ~isempty(obj.cfg.pol)
                fprintf('Polarization : %s\n', obj.cfg.pol);
            else
                fprintf('Polarization : Not set\n');
            end

            % --- Port / Source ---
            if isfield(obj.cfg,'port') && ~isempty(obj.cfg.port)
                p = obj.cfg.port;
                fprintf('Source Port  : dir=%s\n', p.dir);
                % Port1D is on xy
                if p.dir(2) == 'x'
                    fprintf('              x = %.3e m, y-span = [%.3e, %.3e] m\n', p.p1(1), p.p1(2), p.p2(2));
                elseif p.dir(2) == 'y'
                    fprintf('              y = %.3e m, x-span = [%.3e, %.3e] m\n', p.p1(2), p.p1(1), p.p2(1));
                end
            else
                fprintf('Source Port  : Not set\n');
            end

            if isfield(obj.cfg,'source') && ~isempty(obj.cfg.source)
                s = obj.cfg.source;
                if isprop(s,'setflag') && s.setflag
                    fprintf('Source       : Set (type: Source1D, pol: %s)\n', string(s.pol));
                    if isprop(s,'neff') && ~isempty(s.neff)
                        fprintf('              src neff = %s\n', mat2str(s.neff(:).', 6));
                    end
                else
                    fprintf('Source       : Present but not configured\n');
                end
            else
                fprintf('Source       : Not set\n');
            end

            % --- Solution status ---
            if isempty(obj.solution)
                fprintf('Solution     : Not created\n');
            else
                if isprop(obj.solution,'solflag') && obj.solution.solflag
                    fprintf('Solved       : YES\n');
                else
                    fprintf('Solved       : No\n');
                end
            end

            fprintf('==============================\n');
        end

        function results = exportResults(obj)
            if obj.stage < 8
                dispError('Modelvar3D:NoSolutionData');
            end

            results = struct();

            % --- Base metadata ---
            results.wavelength  = obj.lam;
            results.background  = obj.bgmat;
            results.solver      = obj.solver;

            % --- Coordinates ---
                results.x = obj.mesh.axisx.v;
                results.y = obj.mesh.axisy.v;
                results.z = obj.mesh.axisz.v;

            % --- Config  ---
            results.bc = obj.cfg.bc; 
            results.pol = obj.cfg.pol;  
            results.port = obj.cfg.port; 
            results.source = obj.cfg.source;

            % --- varFDFD specific ---
            if isprop(obj.solution,'epseff') && ~isempty(obj.solution.epseff)
                results.epseff = obj.solution.epseff;
            end
            if isprop(obj.solution,'ref') && ~isempty(obj.solution.ref)
                results.ref = obj.solution.ref;  % includes pin, refDevice, mode (EigenMode1D)
            end

            % --- Fields (3D) ---
            results.Ex = obj.solution.Ex;
            results.Ey = obj.solution.Ey;
            results.Ez = obj.solution.Ez;
            results.Hx = obj.solution.Hx;
            results.Hy = obj.solution.Hy;
            results.Hz = obj.solution.Hz;
        end
    end

    methods (Access = private)
        function initvarFDFD(obj, varargin)
            arguments
                obj
            end
            arguments(Repeating)
                varargin
            end

            % parse input
            p = inputParser;
            addParameter(p, 'Pol', 'TE', @(x) ismember(x, {'TE','TM'}));
            addParameter(p, 'RefPoint', [0 0], @(x) isnumeric(x) && numel(x)==2);
            parse(p, varargin{:});

            pol  = p.Results.Pol;
            rpin = p.Results.RefPoint;

            % initialize the var FDFD solver object
            obj.solution = varFDFD(obj.mesh, obj.lam, pol);
            obj.solution = obj.solution.setReference(obj.device, rpin);
        end
    end
end