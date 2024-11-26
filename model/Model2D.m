classdef Model2D < handle
    % Model2D: Class for setting up and solving 2D electromagnetic simulations
    %
    % Description:
    %   The `Model2D` class allows users to set up a two-dimensional electromagnetic model with defined devices,
    %   mesh, ports, sources, and solvers. This class handles the entire simulation workflow, from defining materials and mesh
    %   to adding ports and sources and solving for the electromagnetic fields using FDFD techniques.
    %
    % Properties:
    %   stage - The current stage of the model's configuration.
    %   label - The plane label indicating the propagation plane ('xy', 'yz', 'zx').
    %   lam - The wavelength of the electromagnetic wave in the simulation.
    %   devicelist - List of devices added to the model.
    %   device - The combined device used in the simulation.
    %   mesh - The mesh grid for the simulation.
    %   solver - The type of solver used ('EigenMode' or 'Scattering').
    %   tpml - Thickness of the perfectly matched layers (PML) for the solver.
    %   port - The port used for excitation.
    %   source - The source used for excitation in the model.
    %   solfield - The computed electromagnetic fields after solving the model.
    %
    % Methods:
    %   setLabel(plane_label) - Set the plane label for the simulation.
    %   setWavelength(wavelength) - Set the operating wavelength.
    %   addDevice(device_layer) - Add a device layer to the model.
    %   setDevice(device_layer) - Set the current device for modifications.
    %   finishDevice() - Combine all device layers into a single device.
    %   setMesh(num_point, step_size, varargin) - Set the mesh for the simulation.
    %   setSolver(solver, pml_thickness) - Set the solver type and PML thickness.
    %   addPort(direction, position, span_range, unit) - Add a port to the model.
    %   addSource(source_type, varargin) - Add a source to the model.
    %   setSource() - Get the current source configuration.
    %   solve(solver_parameter) - Solve the model using the configured solver.
    %
    % Example:
    %   % Create an instance of Model2D and set up a simulation
    %   model = Model2D.initialize();
    %   model.setLabel('xy');
    %   model.setWavelength(1.55e-6);
    %   model.addDevice(1);
    %   model.assemblehDevice();
    %   model.setMesh([100, 100], [0.01, 0.01]);
    %   model.setSolver('Scattering', [10, 10]);
    %   model.addPort('+x', 0.5, [0.1, 0.1]);
    %   model.addSource('PlaneWave');
    %   model.solve();
    %
    % Notes:
    %   - This class is intended to be used as a singleton using the `initialize()` method.
    %   - The configuration stages must be followed in order for successful simulation setup.
    %
    % See Also:
    %   Device2D, Source1D, Grid2D, Port1D

    properties (SetAccess = private)
        stage
        label
        lam
        devicelist
        device
        mesh
        solver
        tpml
        port
        source
        solfield
    end

    methods (Access = private)
        function obj = Model2D()
            obj.resetProperties();
        end
    end

    methods

        % Label
        function setLabel(obj, plane_label)
            arguments
                obj
                plane_label {mustBeMember(plane_label, {'xy', 'yz', 'zx'})}
            end

            obj.label = plane_label;
            obj.checkStage([0, 1, 2], 'Model2D:SetLabel');
            obj.updateStageLabelLam;
        end

        % Wavelength
        function setWavelength(obj, wavelength)
            arguments
                obj
                wavelength {mustBePositive}
            end

            obj.lam = wavelength;
            obj.checkStage([0, 1, 2], 'Model2D:SetWavelength');
            obj.updateStageLabelLam();
        end

        % Device
        function addDevice(obj, device_layer)
            arguments
                obj
                device_layer {mustBeInteger, mustBePositive}
            end

            obj.checkStage(2, 'Model2D:AddDevice');
            obj.devicelist(device_layer) = Device2D(device_layer, obj.label);
        end

        function device = setDevice(obj, device_layer)
            arguments
                obj
                device_layer {mustBeInteger, mustBePositive}
            end
            device = obj.devicelist(device_layer);
        end

        function assembleDevice(obj)
            % assembleDevice: Combine all device layers into a single device
            %
            %   Syntax:
            %     finishDevice()
            %
            %   Description:
            %     Combines all device layers into a single device for the simulation.

            combined_device = obj.devicelist(1).copyDevice;
            combined_device = combined_device.combineDevice(obj.devicelist(2:end));
            obj.device = combined_device;
            obj.stage = 3;
        end

        % Mesh
        function setMesh(obj, num_point, step_size, varargin)
            arguments
                obj
                num_point (1,2) double {mustBePositive, mustBeInteger}
                step_size (1,2) double {mustBeReal, mustBePositive}
            end

            arguments(Repeating)
                varargin
            end

            obj.checkStage(3, 'Model2D:SetMesh');
            obj.mesh = Grid2D(num_point, step_size, varargin{:});
            obj.device = obj.device.meshDevice(obj.mesh);
            obj.stage = 4;
        end

        % Solver
        function setSolver(obj, solver, pml_thickness)
            arguments
                obj
                solver {mustBeMember(solver, {'EigenMode', 'Scattering'})}
                pml_thickness = []
            end

            obj.checkStage(4, 'Model2D:SetSolver');
            obj.solver = solver;
            obj.tpml = pml_thickness;
            obj.stage = 5;
        end

        % Port
        function addPort(obj, direction, position, span_range, unit)
            arguments
                obj
                direction {mustBeMember(direction, {'+x', '-x', '+y', '-y', '+z', '-z'})}
                position {mustBeReal}
                span_range (1,2) {mustBeReal}
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            obj.checkStage(5, 'Model2D:AddPort');
            if strcmp(obj.solver, 'EigenMode')
                dispError('Model2D:NoPortNeed');
            end
            if ~contains(obj.label, direction(2))
                dispError('Port1D:DirectionWrong', obj.label(1), obj.label(2));
            end
            obj.port = Port1D(obj.label, direction, position, span_range, unit);
            obj.stage = 6;
        end

        % Source
        function addSource(obj, source_type, varargin)
            arguments
                obj
                source_type {mustBeMember(source_type, {'EigenMode', 'GaussianBeam', 'PlaneWave'})}
            end
            arguments (Repeating)
                varargin
            end

            obj.checkStage(6, 'Model2D:AddSource');
            if strcmp(obj.solver, 'EigenMode')
                dispError('Model2D:NoSourceNeed');
            end

            obj.source = Source1D(obj.label, obj.port, obj.mesh);

            switch source_type
                case 'EigenMode'
                    obj.source = obj.source.setEigenMode(obj.device, obj.lam, varargin{:});
                case 'GaussianBeam'
                    obj.source = obj.source.setGaussianBeam(obj.lam, obj.devicelist(1).mat, varargin{:});
                case 'PlaneWave'
                    obj.source = obj.source.setPlaneWave(obj.lam, obj.devicelist(1).mat, varargin{:});
            end

            obj.stage = 7;
        end

        function source = setSource(obj)
            source = obj.source;
        end

        % FDFD simulation
        function solve(obj, solver_parameter)
            arguments
                obj 
                solver_parameter = []
            end
            if strcmp(obj.solver, 'Scattering') && obj.stage == 7
                solve_result = obj.device.solveScattering(obj.source, obj.tpml);
                obj.solfield.Ex = solve_result.Ex;
                obj.solfield.Ey = solve_result.Ey;
                obj.solfield.Ez = solve_result.Ez;
                obj.solfield.Hx = solve_result.Hx;
                obj.solfield.Hy = solve_result.Hy;
                obj.solfield.Hz = solve_result.Hz;
                obj.stage = 8;
            elseif strcmp(obj.solver, 'EigenMode') && obj.stage == 5
                % 'solver_parameter' is mode number
                solve_result = obj.device.solveEigenMode(obj.lam, solver_parameter);
                obj.solfield.Ex = solve_result.Ex;
                obj.solfield.Ey = solve_result.Ey;
                obj.solfield.Ez = solve_result.Ez;
                obj.solfield.Hx = solve_result.Hx;
                obj.solfield.Hy = solve_result.Hy;
                obj.solfield.Hz = solve_result.Hz;
                obj.stage = 6;
            else
                dispError('Model2D:Solve');
            end
        end

        function results = exportResult(obj)
            switch obj.solver
                case 'Scattering'
                    checkStage(obj, 8, 'Model2D:ScatteringNotSolved');
                case 'EigenMode'
                    checkStage(obj, 6, 'Model2D:EigenModeNotSolved');
            end

            results.(sprintf(obj.label(1))) = obj.mesh.axis1.v;
            results.(sprintf(obj.label(2))) = obj.mesh.axis2.v;
            results.Ex = obj.solfield.Ex;
            results.Ey = obj.solfield.Ey;
            results.Ez = obj.solfield.Ez;
            results.Hx = obj.solfield.Hx;
            results.Hy = obj.solfield.Hy;
            results.Hz = obj.solfield.Hz;
        end

        function dispField(obj, field_name, data_type)
            arguments
                obj 
                field_name {mustBeMember(field_name, {'Ex','Ey','Ez','Hx','Hy','Hz'})}
                data_type {mustBeMember(data_type, {'abs', 'real', 'imag'})} = 'abs'
            end

            field = obj.solfield.(field_name);
            switch data_type
                case 'abs'
                    field_data = abs(field);
                case 'real'
                    field_data = real(field);
                case 'imag'
                    field_data = imag(field);
            end

            for ii = 1:size(field,3)
                figure;
                imagesc(obj.mesh.axis1.v, obj.mesh.axis2.v, field_data(:,:,ii)');
                xlabel(obj.label(1)), ylabel(obj.label(2)), axis image;
                title(['Field ',num2str(ii), ' , ', data_type, '(', field_name, ')']);
            end
        end
    end

    methods (Static)
        function singleObj = initialize()
            persistent uniqueInstance

            if isempty(uniqueInstance) || ~isvalid(uniqueInstance)
                uniqueInstance = Model2D();
                disp("Model2D: Object Initialized"); % Display message for a new instance
            else
                % Reset properties of the existing instance
                uniqueInstance.resetProperties();
                disp("Model2D: Object Reset and Initialized"); % Display message for resetting an existing instance
            end

            singleObj = uniqueInstance;
        end
    end

    methods (Access = private)

        function resetProperties(obj)
            % Reset all properties to their initial values
            obj.stage = 0;
            obj.label = [];
            obj.lam = [];
            obj.devicelist = Device2D.empty(0,1);
            obj.mesh = [];
            obj.solver = [];
            obj.solfield = [];
        end

        function checkStage(obj, allowedStages, error_id)
            % Check whether the current stage is one of the allowed stages
            if ~ismember(obj.stage, allowedStages)
                dispError(error_id);
            end
        end

        function updateStageLabelLam(obj)
            % Update the stage for 'label' and 'lam' properties.
            if obj.stage == 0
                obj.stage = 1;  % First property (either label or lam) set.
            elseif obj.stage == 1
                obj.stage = 2;  % Second property (either label or lam) set.
            end
        end

    end
end