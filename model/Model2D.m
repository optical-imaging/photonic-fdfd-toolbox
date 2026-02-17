classdef Model2D < Model
    %Model2D: 2D model workflow for building a Device2D + Grid2D and running Scattering2D or EigenMode2D
    %
    % Key Properties (inherited from Model, SetAccess = protected):
    %   stage       - Workflow stage/state
    %   lam         - Wavelength (m)
    %   bgmat       - Background material (Material)
    %   devicelist  - Device2D component list (indexed by seq)
    %   device      - Assembled Device2D used for simulation
    %   mesh        - Grid2D mesh
    %   solver      - Solver selection: 'Scattering'|'EigenMode'
    %   cfg         - Configuration struct (pol, bc, port, source, num_mode, etc.)
    %   solution    - Solver object (Scattering2D or EigenMode2D) holding solved fields
    %
    % Key Properties (Model2D-specific):
    %   label       - Simulation plane label: 'xy'|'yz'|'zx' (defines axis ordering)
    %
    % Key Methods (Initialization / reset):
    %   Model2D()                       - Constructor (called only by Model.initializeSingleObj)
    %   resetProperties()               - Reset Model2D state and allocate empty Device2D containers
    %   initialize()                    - Return singleton Model2D instance
    %
    % Key Methods (Geometry / mesh):
    %   setLabel(plane_label)           - Set simulation plane label before meshing
    %   setMesh(num_point,step_size)    - Create Grid2D mesh (overrides Model.setMesh signature to 2D)
    %   getCoordinate()                 - Return axis coordinate vectors from current Grid2D mesh
    %   addDeviceBit(eps_bit,mu_bit,step) - Define Device2D directly from 2D bitmap arrays
    %
    % Key Methods (Solver selection and configuration):
    %   setSolver(solver_type,var3)     - Create and configure solution object:
    %                                    - 'Scattering': var3 is polarization 'TE'|'TM', default BC is PML/PML
    %                                    - 'EigenMode':  var3 is number of modes (default 1), default BC is Perfect/Perfect
    %   setBoundary(axis_label,type,opt)- Update BC type for one axis; for PML, opt sets thickness
    %
    % Key Methods (Ports and sources; Scattering-only):
    %   addPort(direction,position,span)- Create Port1D on the current plane (calls Model.addPort)
    %   addSource(source_type,opt)      - Build Source1D on cfg.port and attach to Scattering2D:
    %                                    - 'EigenMode': opt is mode indices/count passed to Source1D.setEigenMode
    %                                    - 'GaussianBeam': opt is waist
    %                                    - 'PlaneWave': no opt
    %   setSource()                     - Return cfg.source handle (inherited helper)
    %   getSource()                     - Return cfg.source if exists (inherited helper)
    %
    % Key Methods (Run and output):
    %   solveModel()                    - Execute selected solver:
    %                                    - Scattering: requires port+source stage, calls Scattering2D.solveFDFD(device)
    %                                    - EigenMode: calls EigenMode2D.solveFDFD(device)
    %   dispInfo()                      - Print current configuration summary (developer debugging)
    %   exportResults()                 - Export coordinates + fields, plus bc (Scattering) or neff/mode_index (EigenMode)

    properties (SetAccess = private)
        label
    end

    methods (Access = private)
        function obj = Model2D()
            obj.resetProperties();
        end

        function resetProperties(obj)
            obj.resetPartialProperties; % inherite from Model
            obj.label = 'xy';
            obj.devicelist = Device2D.empty(0,1);
            obj.device = Device2D.empty(0,1);
        end
    end

    methods (Static)
        function singleObj = initialize()
            persistent uniqueInstance

            if isempty(uniqueInstance) || ~isvalid(uniqueInstance)
                disp("Model2D: Object Initialized");
            else
                disp("Model2D: Object Reset and Initialized");
            end

            uniqueInstance = Model2D();
            singleObj = uniqueInstance;
        end
    end

    methods
        % Label
        function setLabel(obj, plane_label) % optional
            arguments
                obj
                plane_label {mustBeMember(plane_label, {'xy', 'yz', 'zx'})}
            end

            obj.checkStage([0, 1], 'Model2D:ChangeLabel');
            obj.label = plane_label;
        end

        % Mesh
        function setMesh(obj, num_point, step_size)
            arguments
                obj
                num_point (1,2) double {mustBePositive, mustBeInteger}
                step_size (1,2) double {mustBeReal, mustBePositive}
            end
            obj.checkStage([2,3,4], 'Model:SetMesh');

            obj.mesh = Grid2D(num_point, step_size, obj.label);
            obj.stage = 3;
        end

        function [axis1_v, axis2_v] = getCoordinate(obj)
            if isempty(obj.mesh)
                dispError('Model:NoMesh');
            end
            axis1_v = obj.mesh.axis1.v;
            axis2_v = obj.mesh.axis2.v;
        end

        % Directly set bitmap
        function addDeviceBit(obj, eps_bit, mu_bit, step_size)
            arguments
                obj
                eps_bit (:,:) double {mustBeRealAndGE1}
                mu_bit (:,:) double {mustBeRealAndGE1}
                step_size (1,2) double {mustBeReal, mustBePositive}
            end

            obj.checkStage(1, 'Model:AddDevice');

            obj.device = Device2D(1);
            obj.device.setBitmap(eps_bit, mu_bit, step_size, obj.label); % 2D grid needs label

            obj.mesh = obj.device.mesh;

            obj.stage = 4;
        end

        % Solver
        function setSolver(obj, solver_type, var3)
            arguments
                obj
                solver_type {mustBeMember(solver_type, {'EigenMode', 'Scattering'})}
                var3 = [] % Default to empty
            end

            obj.checkStage([4,5,6,7,8], 'Model2D:SetSolver');

            if isempty(obj.lam)
                dispError('Model:SetWavelength');
            end

            obj.solver = solver_type;
            obj.cfg = struct(); % reset cfg

            switch solver_type
                case 'Scattering'
                    % Validation: Must be TE or TM
                    if ~ismember(var3, {'TE', 'TM'})
                        dispError('Model2D:SetScatteringPol');
                    end

                    obj.cfg.pol = var3; % Store polarization

                    % PML is default boundary condition for Scattering solver
                    bc1 = BC1D(obj.label(1), 'PML'); % default PML has 8 layers
                    bc2 = BC1D(obj.label(2), 'PML');
                    obj.cfg.bc = BC2D(bc1, bc2);

                    % set up Scattering2D solver
                    obj.solution = Scattering2D(obj.lam, obj.cfg.pol);
                    obj.solution.setMesh(obj.mesh);
                    obj.solution.setBC(obj.cfg.bc);

                    obj.stage = 5;

                case 'EigenMode'
                    % Validation: Must be a positive integer
                    if isempty(var3)
                        num_mode = 1; % Default value
                    elseif isnumeric(var3) && isscalar(var3) && var3 > 0 && mod(var3, 1) == 0
                        num_mode = var3;
                    else
                        error('For EigenMode, the 3rd argument must be a positive integer.');
                    end

                    obj.cfg.num_mode = num_mode;

                    % Perfect boundary condition is default for EigenMode2D
                    bc1 = BC1D(obj.label(1), 'Perfect');
                    bc2 = BC1D(obj.label(2), 'Perfect');
                    obj.cfg.bc = BC2D(bc1, bc2);

                    % set up Scattering2D solver
                    obj.solution = EigenMode2D(obj.lam, num_mode);
                    obj.solution.setMesh(obj.mesh);
                    obj.solution.setBC(obj.cfg.bc);

                    obj.stage = 7;
            end

        end

        % Boundary condition
        function setBoundary(obj, axis_label, BC_type, opt)
            arguments
                obj
                axis_label {mustBeMember(axis_label, {'x','y','z'})}
                BC_type {mustBeMember(BC_type, {'Perfect', 'PML', 'Periodic'})}
                opt = [] % for PML to set number of layers
            end

            obj.checkStage([5,6,7,8], 'Model2D:SetBC');

            % check axis name
            if ~contains(obj.label, axis_label)
                dispError('Model2D:BCaxisWrong', obj.label(1), obj.label(2));
            end

            % Solver-specific policy
            if strcmp(obj.solver, 'EigenMode') && strcmp(BC_type, "PML")
                dispError('Model:EigenModeNoPML');   % can't assign PML for EigenMode solver
            end
            if strcmp(obj.solver, 'Scattering') && strcmp(BC_type, "Perfect")
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

        % Port
        function addPort(obj, direction, position, span_range)
            arguments
                obj
                direction {mustBeMember(direction, {'+x', '-x', '+y', '-y', '+z', '-z'})}
                position double {mustBeReal}
                span_range (1,2) double {mustBeReal}
            end

            % Eigenmode solver doesn't need port
            if strcmp(obj.solver, 'EigenMode')
                dispError('Model2D:NoPortNeed');
            end

            obj.checkStage([5, 6, 7, 8], 'Model:AddPort');
            
            % Port has be in plane
            if ~contains(obj.label, direction(2))
                dispError('Port1D:DirectionWrong', obj.label(1), obj.label(2));
            end

            obj.cfg.port = Port1D(obj.label, direction, position, span_range);
            obj.cfg.source = []; % clear source if port has been changed

            obj.stage = 6;
        end

        % Source
        function addSource(obj, source_type, opt)
            arguments
                obj
                source_type {mustBeMember(source_type, {'EigenMode', 'GaussianBeam', 'PlaneWave'})}
                opt = []
            end

            if strcmp(obj.solver, 'EigenMode')
                dispError('Model2D:NoSourceNeed');
            end

            obj.checkStage([6,7,8], 'Model2D:AddSource');

            s = Source1D(obj.lam, obj.cfg.port, obj.cfg.pol);
            s.setMesh(obj.mesh);
            switch source_type
                case 'EigenMode'
                    s.setEigenMode(obj.device, opt); % opt is number of modes for EigenMode solver
                case 'GaussianBeam'
                    s.setGaussianBeam(obj.bgmat, opt); % opt is waist width for GaussianBeam
                case 'PlaneWave'
                    s.setPlaneWave(obj.bgmat);
            end

            obj.cfg.source = s;
            obj.solution.setSrc(s);

            obj.stage = 7;
        end

        % solve
        function solveModel(obj)
            obj.checkStage(7, 'Model2D:SetBeforeSolving');

            switch obj.solver
                case 'Scattering'
                    obj.solution.solveFDFD(obj.device);
                case 'EigenMode'
                    obj.solution.solveFDFD(obj.device);
                otherwise
                    dispError('Model2D:NoSolverSelected');
            end

            obj.stage = 8;
            dispMessage('Model:ModelSolved');
        end

        % export
        function dispInfo(obj)
            % dispInfo  Show a summary of the current Model2D configuration
            fprintf('=== Model2D Information ===\n');

            % --- General ---
            fprintf('Wavelength   : %.3e m\n', obj.lam); % Using scientific notation
            fprintf('Plane Label  : %s\n', obj.label);
            fprintf('Background   : eps_r=%.3f, mu_r=%.3f\n', obj.bgmat.eps, obj.bgmat.mu);

            % --- Mesh ---
            if isempty(obj.mesh)
                fprintf('Mesh         : Not set\n');
            else
                ax1 = obj.mesh.axis1;
                ax2 = obj.mesh.axis2;
                lbl1 = obj.label(1);
                lbl2 = obj.label(2);
                % Removed .u property and changed step size to scientific notation %.3e
                fprintf('Mesh         : %s: N=%d, d=%.3e m\n', lbl1, ax1.n, ax1.d);
                fprintf('               %s: N=%d, d=%.3e m\n', lbl2, ax2.n, ax2.d);
            end

            % --- Solver ---
            if isempty(obj.solver)
                fprintf('Solver       : Not set\n');
            else
                fprintf('Solver       : %s\n', obj.solver);
            end

            % --- Solver-specific ---
            switch obj.solver
                case 'Scattering'
                    % Boundary conditions
                    if isfield(obj.cfg,'bc') && ~isempty(obj.cfg.bc)
                        bc = obj.cfg.bc;
                        for b = {'bc1','bc2'}
                            if isprop(bc,b{1})
                                bci = bc.(b{1});
                                if strcmp(bci.type,'PML')
                                    fprintf('BC (%s)     : %s, thickness=%d\n', ...
                                        bci.label, bci.type, bci.pmlnum);
                                else
                                    fprintf('BC (%s)     : %s\n', bci.label, bci.type);
                                end
                            end
                        end
                    else
                        fprintf('Boundary Cond: Not set\n');
                    end

                    % Polarization
                    if isfield(obj.cfg,'source') && ~isempty(obj.cfg.source)
                        fprintf('Source Polarization : %s\n', obj.cfg.source.pol);
                    else
                        fprintf('Source Polarization : Not set\n');
                    end

                    % Port
                    if isfield(obj.cfg,'port') && ~isempty(obj.cfg.port)
                        p = obj.cfg.port;
                        switch p.dir(2)
                            case obj.label(1)
                                pos_axis = obj.label(1);
                                span_axis = obj.label(2);
                                position = p.p1(1);
                                span = [p.p1(2), p.p2(2)];
                            case obj.label(2)
                                pos_axis = obj.label(2);
                                span_axis = obj.label(1);
                                position = p.p1(2);
                                span = [p.p1(1), p.p2(1)];
                            otherwise
                                pos_axis = '?'; span_axis = '?';
                                position = NaN; span = [NaN NaN];
                        end
                        fprintf('Source Port         : dir=%s\n', p.dir);
                        % Switched position and span to %.4e for better resolution in meters
                        fprintf('                      position (%s) = %.3e m\n', pos_axis, position);
                        fprintf('                      span (%s)     = [%.3e, %.3e] m\n', span_axis, span(1), span(2));
                    else
                        fprintf('Source Port         : Not set\n');
                    end

                case 'EigenMode'
                    if isfield(obj.cfg,'num_mode') && ~isempty(obj.cfg.num_mode)
                        fprintf('Requested Modes: %s\n', mat2str(obj.cfg.num_mode));
                    end
                    if ~isempty(obj.solution) && isprop(obj.solution,'neff')
                        fprintf('Modes and neff:\n');
                        for k = 1:numel(obj.solution.neff)
                            neff_val = obj.solution.neff(k);
                            if abs(imag(neff_val)) < 1e-6
                                fprintf('   Mode %d : neff = %.6f\n', obj.cfg.num_mode(k), real(neff_val));
                            else
                                fprintf('   Mode %d : neff = %.6f %+.6fi\n', ...
                                    obj.cfg.num_mode(k), real(neff_val), imag(neff_val));
                            end
                        end
                    end
            end
            fprintf('===========================\n');
        end

        function results = exportResults(obj)
            if isempty(obj.solution)
                dispError('Model2D:NoSolutionData');
            end

            % --- Coordinates ---
            axis1 = obj.mesh.axis1.v;
            axis2 = obj.mesh.axis2.v;

            % --- Base metadata ---
            results = struct();
            results.wavelength    = obj.lam;
            results.background  = obj.bgmat;
            results.solver = obj.solver;

            % --- Assign coordinates directly by label ---
            results.(obj.label(1)) = axis1;
            results.(obj.label(2)) = axis2;

            % --- Assign all field components directly ---
            results.Ex = obj.solution.Ex;
            results.Ey = obj.solution.Ey;
            results.Ez = obj.solution.Ez;
            results.Hx = obj.solution.Hx;
            results.Hy = obj.solution.Hy;
            results.Hz = obj.solution.Hz;

            % --- Solver-specific additions ---
            switch obj.solver
                case 'Scattering'
                    results.bc = obj.cfg.bc;

                case 'EigenMode'
                    results.neff = obj.solution.neff;
                    results.mode_index = obj.cfg.num_mode;
            end
        end
    end

end