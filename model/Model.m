classdef (Abstract, Hidden) Model < handle
    %Model: Abstract workflow controller that builds Device/Mesh/BC/Source and runs an FDFD solver
    %
    % Key Properties (SetAccess = protected):
    %   stage       - Workflow stage/state integer controlling allowed operations
    %   lam         - Simulation wavelength (m)
    %   bgmat       - Background material (Material)
    %   devicelist  - List of device components (Device2D or LayerDevice entries indexed by seq)
    %   device      - Assembled simulation device (Device2D | LayerDevice | Device3D)
    %   mesh        - Simulation mesh (Grid2D | Grid3D)
    %   solver      - Solver object (Scattering*/EigenMode* depending on subclass configuration)
    %   cfg         - Configuration struct (ports, sources, BC settings, etc.)
    %   solution    - Solution/results container (structure depends on subclass/export policy)
    %
    % Key Methods (Abstract):
    %   getCoordinate(...)        - Return coordinate arrays/axes (output differs by 2D vs 3D)
    %   setBoundary(...)          - Configure boundary conditions (2D EigenMode path may not require source BC)
    %   addSource(...)            - Create/attach a source (Modelvar3D may require effective index/device)
    %   solveModel(...)           - Run solver workflow (device meshing, solver setup, solve, postprocess)
    %   dispInfo(...)             - Display model/solver status (varies by subclass)
    %   exportResults(...)        - Export results package (content varies by subclass)
    %
    % Key Methods (Global settings):
    %   setWavelength(lam)        - Set wavelength and advance stage from 0->1 if needed
    %   setBackground(opts)       - Set background material via Material or (eps_r,mu_r)
    %
    % Key Methods (Device construction):
    %   addDevice(seq)            - Create/append a device component at index seq (enforces contiguous 1..K)
    %   setDevice(seq)            - Return handle to device component for configuration
    %   assembleDevice()          - Combine devicelist into a single simulation device and set stage=2
    %
    % Key Methods (Meshing):
    %   setMesh(num_point,step)   - Create Grid2D/Grid3D mesh and set stage=3
    %   meshDevice()              - Rasterize assembled device onto current mesh and set stage=3
    %   shiftMesh(dim,mode,val)   - Shift mesh grid coordinates (delegates to mesh.shiftGrid)
    %
    % Key Methods (Bitmap-defined device):
    %   addDeviceBit(eps_bit,mu_bit,step) - Directly define device from bitmap arrays and set stage=3
    %
    % Key Methods (Port/Source access):
    %   addPort(dir,pos,span)     - Create Port1D/Port2D and clear existing source, set stage=6
    %   setSource()               - Return current cfg.source handle for configuration
    %   getSource()               - Return cfg.source if exists; otherwise []
    %
    % Key Methods (Static, Access = protected):
    %   initializeSingleObj(ndim) - Return singleton Model instance for '2D'|'Semi3D'|'3D' and reset if reused
    %
    % Key Methods (Access = protected):
    %   resetPartialProperties()  - Reset core workflow properties to initial defaults
    %   checkStage(stages,id)     - Enforce allowed stage transitions (throws via dispError)
    %   checkDeviceIndex(idx)     - Validate that devicelist(idx) exists and is non-empty
    %   resetDevice()             - Clear assembled device/mesh/solver/cfg/solution and return to stage=1

    properties (SetAccess = protected)
        stage
        lam
        bgmat
        devicelist
        device
        mesh
        solver
        cfg
        solution
    end

    methods (Abstract)
        setMesh
        addDeviceBit
        getCoordinate % output is different for 2D and 3D'
        setBoundary % Model2D Eigenmode solver doesn't need
        addPort
        addSource % Modelvar3D need effective index as device
        solveModel % solving configurations varies
        dispInfo % solvers info varies
        exportResults % different data will be included
    end

    methods
        %% global setting
        % Wavelength
        function setWavelength(obj, wavelength) % required
            arguments
                obj
                wavelength double {mustBePositive} % unit: m
            end

            obj.lam = wavelength;

            if obj.stage >= 5 % after solver has been set
                dispWarning('Model:WavelengthChange');
                obj.stage = 4; % just before solver configuration
            end
        end

        % Background
        function setBackground(obj, varargin)
            arguments
                obj
            end
            arguments (Repeating)
                varargin
            end

            obj.checkStage([1,2,3], 'Model:ChangeBackground');

            % input Material
            if isscalar(varargin) && isa(varargin{1}, 'Material')
                obj.bgmat = varargin{1};
                return;
            end

            % input values
            p = inputParser;
            p.FunctionName = 'Model.setBackground';

            addParameter(p, 'eps_r', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(p, 'mu_r',  1, @(x) isnumeric(x) && isscalar(x) && x > 0);

            parse(p, varargin{:});

            eps_r = p.Results.eps_r;
            mu_r  = p.Results.mu_r;

            obj.bgmat = Material(eps_r, mu_r);
        end

        %% sumulation device
        % initialize device components
        function addDevice(obj, device_seq)
            arguments
                obj
                device_seq {mustBeInteger, mustBePositive}
            end

            obj.checkStage([1,2,3,4], 'Model:AddDevice');
            if obj.stage >= 2
                obj.resetDevice();
            end

            % Existing present indices (contiguous 1..K by construction below)
            K = 0;
            if ~isempty(obj.devicelist)
                present_mask = false(1, numel(obj.devicelist));
                for k = 1:numel(obj.devicelist)
                    present_mask(k) = ~isempty(obj.devicelist(k));
                end
                K = sum(present_mask);  % since we enforce contiguity, present_idx = 1:K
            end

            submodel = class(obj);
            if obj.checkDeviceIndex(device_seq)
                switch submodel
                    case 'Model2D'
                        obj.devicelist(device_seq) = Device2D(device_seq); % Replace/update existing
                    case 'Modelvar3D'
                        obj.devicelist(device_seq) = LayerDevice(device_seq);
                end
            else
                if device_seq ~= K + 1 % Append new layer — must be exactly next index K+1 (no gaps)
                    dispError('Model:DeviceIndexGap');  % must be 1,2,3,...,K,(K+1)
                end
                switch submodel % Expand array safely to device_seq and assign
                    case 'Model2D'
                        obj.devicelist(device_seq) = Device2D(device_seq); % Replace/update existing
                    case 'Modelvar3D'
                        obj.devicelist(device_seq) = LayerDevice(device_seq);
                end
            end
        end

        % config single device component
        function device = setDevice(obj, device_seq)
            arguments
                obj
                device_seq {mustBeInteger, mustBePositive}
            end

            obj.checkStage([1,2,3], 'Model:AddDevice');
            if obj.stage >= 2
                obj.resetDevice();
            end

            if ~obj.checkDeviceIndex(device_seq)
                dispError('Model:DeviceMissing');  % cannot select non-existent/gap
            end

            device = obj.devicelist(device_seq); % handle pass to be set
        end

        % combine device components as a single device for simulation
        function assembleDevice(obj)
            if isempty(obj.devicelist)
                dispError('Model:NoDevice');
            end

            empties     = arrayfun(@isempty, obj.devicelist);
            present_idx = find(~empties);
            if isempty(present_idx) || ~isequal(present_idx, 1:present_idx(end))
                dispError('Model:DeviceIndexGap');  % Either all entries are empty, or indices are not 1..k
            end

            if numel(obj.devicelist) > 1
                if ~isa(obj, 'Modelvar3D')
                    combined_device = Device.combineDevice(obj.devicelist(:));
                else
                    combined_device = LayerDevice.combineDevice(obj.devicelist(:));
                end
            else
                combined_device = obj.devicelist(1);
            end

            obj.device = combined_device;
            obj.stage = 2;
        end

        %% device meshing
        % mesh Device object as bit array
        function meshDevice(obj)
            obj.checkStage(3, 'Model:MeshDevice')
            if isempty(obj.mesh)
                dispError('Model:NoMesh');
            end
            
            obj.device.meshDevice(obj.mesh, obj.bgmat);
            obj.stage = 4;
        end

        % customize mesh grid
        function shiftMesh(obj, dim, mode, new_value)
            if isempty(obj.mesh)
                dispError('Model:NoMesh');
            end
            obj.mesh = obj.mesh.shiftGrid(dim, mode, new_value);
        end

        %% port & source
        % config source
        function source = setSource(obj)
            source = obj.cfg.source;
        end

        % get source object (if existing)
        function source = getSource(obj)
            if isfield(obj.cfg,'source')
                source = obj.cfg.source;
            else
                source = [];
            end
        end

    end

    methods (Access = protected)
        function resetPartialProperties(obj)
            obj.stage = 1;
            obj.lam = [];
            obj.bgmat    = Material(1,1);
            obj.mesh = [];
            obj.solver = [];
            obj.cfg = struct();
            obj.solution = [];
        end

        function checkStage(obj, allowedStages, error_id)
            if ~ismember(obj.stage, allowedStages)
                dispError(error_id);
            end
        end

        function tf = checkDeviceIndex(obj, idx)
            tf = (isscalar(idx) && idx >= 1 && idx <= numel(obj.devicelist));
            if tf
                tf = ~isempty(obj.devicelist(idx));
            end
        end

        function resetDevice(obj)
            submodel = class(obj);
            switch submodel
                case 'Model2D'
                    obj.device   = Device2D.empty(0,1);  % assembled device cleared
                case 'Modelvar3D'
                    obj.device   = LayerDevice.empty(0,1);
                case 'Model3D'
                    obj.device   = Device3D.empty(0,1);
            end

            if obj.stage >= 3
                obj.mesh = [];
            end

            obj.stage    = 1; % Back to device set up
            obj.solver   = [];
            obj.cfg      = struct();
            obj.solution = [];
        end

    end

end