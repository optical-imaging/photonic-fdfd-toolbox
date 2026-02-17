classdef Device2D < Device
    %Device2D: 2D simulation device (planar geometry meshed onto Grid2D)
    %
    % Key Properties (inherited, 2D meaning):
    %   mesh       - Grid2D
    %   eps, mu    - 2D maps (expanded internally to include 3 components)
    %
    % Key Methods:
    %   setGeometry(...)        - Assign 2D Layout geometry into geo
    %   meshDevice(grid2d)      - Rasterize geometry onto Grid2D

    methods
        % constructor
        function obj = Device2D(device_seq)
            arguments
                device_seq (1,1) {mustBeInteger, mustBePositive}
            end
            obj@Device(device_seq);
        end

        % device set up
        function setGeometry(obj, shape_name, varargin)
            arguments
                obj
                shape_name {mustBeMember(shape_name, {'Rectangle', 'Disk', 'Polygon', 'Ring', 'GDSII', 'Bitmap'})}
            end
            arguments (Repeating)
                varargin
            end

            switch shape_name
                case 'Rectangle'
                    obj.geo = Rectangle(varargin{:});
                case 'Disk'
                    obj.geo = Disk(varargin{:});
                case 'Polygon'
                    obj.geo = Polygon(varargin{:});
                case 'Ring'
                    obj.geo = Ring(varargin{:});
                case 'GDSII'
                    obj.geo = GDSII(varargin{:});
            end
        end

        % device meshing
        function meshDevice(obj, grid, bg_mat)
            arguments
                obj
                grid Grid2D
                bg_mat Material = []
            end

            if isempty(bg_mat)
                eps_bg = 1;
                mu_bg = 1;
            else
                eps_bg = bg_mat.eps;
                mu_bg = bg_mat.mu;
            end

            % double-sampled grid for Yee cell averaging
            axis1_ds = grid.axis1.doublesampleAxis;
            axis2_ds = grid.axis2.doublesampleAxis;
            [AXIS1, AXIS2] = ndgrid(axis1_ds.v, axis2_ds.v);

            eps_ds = eps_bg*ones(axis1_ds.n, axis2_ds.n);
            mu_ds  = mu_bg*ones(axis1_ds.n, axis2_ds.n);

            % Normalize geometry list
            if isempty(obj.geo)
                geolist = {};
            elseif iscell(obj.geo)
                geolist = obj.geo;
            elseif isvector(obj.geo)
                geolist = num2cell(obj.geo);
            else
                geolist = {obj.geo};
            end

            % Iterate through geometries and assign material values
            for ii = 1:numel(geolist)
                g = geolist{ii};
                
                % find indices of where the geometry be
                in = obj.findInGeometry(AXIS1, AXIS2, g);
                
                eps_val = 1; mu_val  = 1;
                if ii <= numel(obj.mat)
                    m = obj.mat(ii);
                    if isnumeric(m)
                        eps_val = m;
                    else
                        if isprop(m,'eps') && ~isempty(m.eps), eps_val = m.eps; end
                        if isprop(m,'mu')  && ~isempty(m.mu),  mu_val  = m.mu;  end
                    end
                end
                eps_ds(in) = eps_val;
                mu_ds(in)  = mu_val;
            end

            % Sample double-grid maps back to Yee grid components
            [eps1, eps2, eps3, mu1, mu2, mu3] = obj.sampleMapAll(eps_ds, mu_ds);
            obj.eps = cat(3, eps1, eps2, eps3);
            obj.mu  = cat(3, mu1, mu2, mu3);
            obj.mesh = grid;
        end

        function setBitmap(obj, eps_bit, mu_bit, step_size, plane_label)
            arguments
                obj
                eps_bit (:,:) double {mustBeRealAndGE1}
                mu_bit (:,:) double {mustBeRealAndGE1}
                step_size (1,2) double {mustBeReal, mustBePositive}
                plane_label string {mustBeMember(plane_label,{'xy','yz','zx'})}
            end
            % eps and mu size match
            if ~isequal(size(eps_bit), size(mu_bit))
                dispError('Device:BitmapSizeMismatch');
            end

            obj.eps = repmat(eps_bit, [1 1 3]);
            obj.mu  = repmat(mu_bit, [1 1 3]);
            obj.mesh = Grid2D(size(eps_bit), step_size, plane_label);
        end

        % Display
        function dispImg(obj, plot_part, plane_label, varargin)
            % Examples:
            %   obj.dispImg('eps');                                      % auto bbox
            %   obj.dispImg('mu', 'x', [0 5e-6]);                        % set physical x only
            %   obj.dispImg('eps', 'z', [-1e-6 3e-6], 'x', [0 4e-6]);    % label-aware limits
            arguments
                obj
                plot_part {mustBeMember(plot_part, {'eps','mu'})} = 'eps'
                plane_label {mustBeMember(plane_label, {'xy','yz','zx'})} = 'xy'
            end
            arguments (Repeating)
                varargin
            end

            % -------- Parse optional axis limits: reqLims.x / reqLims.y / reqLims.z --------
            reqLims = struct();
            k = 1;
            while k <= numel(varargin)
                key = lower(string(varargin{k}));
                if ~ismember(key, ["x","y","z"])
                    error('dispImg:BadArg', 'Unknown axis "%s". Use ''x'',''y'',''z''.', key);
                end
                if k == numel(varargin) || ~isnumeric(varargin{k+1}) || numel(varargin{k+1}) ~= 2
                    error('dispImg:BadArg', 'Limits for axis %s must be a 1x2 numeric vector.', key);
                end
                reqLims.(key) = varargin{k+1};
                k = k + 2;
            end

            % -------- Decode obj.label and map physical axes -> plotting axes --------
            lbl = lower(plane_label);
            if ~(ischar(lbl) || isstring(lbl)) || numel(lbl) ~= 2 || ~all(ismember(lbl,'xyz'))
                error('dispImg:Label', 'obj.label must be a 2-char string like ''xy'',''yz'',''zx''.');
            end
            axH = lbl(1);  % physical axis shown on plotting X
            axV = lbl(2);  % physical axis shown on plotting Y

            rq = fieldnames(reqLims);
            bad = setdiff(rq, {axH, axV});
            if ~isempty(bad)
                error('dispImg:AxisNotInLabel', ...
                    'Requested axis "%s" not present in label "%s".', bad{1}, obj.label);
            end

            % -------- Normalize obj.geo to a cell array of geometry objects --------
            if isempty(obj.geo)
                geolist = {};
            elseif iscell(obj.geo)
                geolist = obj.geo;
            elseif isvector(obj.geo)
                geolist = num2cell(obj.geo);      % e.g., object array -> cell
            else
                geolist = {obj.geo};              % single object
            end

            % -------- Compute auto bounding box from polygon vertices --------
            if ~isempty(geolist)
                try
                    allx = cellfun(@(g) g.s.Vertices(:,1), geolist, 'UniformOutput', false);
                    ally = cellfun(@(g) g.s.Vertices(:,2), geolist, 'UniformOutput', false);
                catch
                    % If your geometry stores the polyshape differently, adapt here.
                    error('dispImg:NoVertices', ...
                        'Geometry elements must expose polyshape in field ".s".');
                end
                autoX = [min(cellfun(@min, allx)), max(cellfun(@max, allx))];
                autoY = [min(cellfun(@min, ally)), max(cellfun(@max, ally))];
                pad = 0.01 * max([diff(autoX), diff(autoY), eps]); % small margin
                autoX = autoX + [-pad pad];
                autoY = autoY + [-pad pad];
            else
                autoX = [0 1];
                autoY = [0 1];
            end

            % -------- Final plotting ranges (manual overrides win) --------
            xRange = autoX;
            yRange = autoY;
            if isfield(reqLims, axH), xRange = reqLims.(axH); end
            if isfield(reqLims, axV), yRange = reqLims.(axV); end

            % -------- Colormap & figure setup --------
            cmin = 1; cmax = 15; cmap = jet(256);
            figure; hold on;

            % Title
            switch plot_part
                case 'eps', ttl = '\epsilon_r';
                case 'mu',  ttl = '\mu_r';
            end
            in1 = inputname(1);
            if ~isempty(in1), title([in1 ' ' ttl], 'Interpreter','tex');
            else,             title(ttl, 'Interpreter','tex');
            end

            % Background (air = 1)
            imagesc([xRange(1) xRange(2)], [yRange(1) yRange(2)], ones(2));
            set(gca,'YDir','normal');

            % -------- Draw polygons with material colors --------
            N = numel(geolist);
            for ii = 1:N
                % Handle both Material-object and numeric storage
                m = [];
                if ii <= numel(obj.mat), m = obj.mat(ii); end

                switch plot_part
                    case 'eps'
                        if isnumeric(m)
                            val = real(m);
                        elseif ~isempty(m) && isprop(m,'eps') && ~isempty(m.eps)
                            val = real(m.eps);
                        else
                            val = 1;
                        end
                    case 'mu'
                        if isnumeric(m)
                            val = real(m);
                        elseif ~isempty(m) && isprop(m,'mu') && ~isempty(m.mu)
                            val = real(m.mu);
                        else
                            val = 1;
                        end
                end

                t = (val - cmin) / (cmax - cmin);
                t = max(0, min(1, t));
                faceCol = cmap(round(t*(size(cmap,1)-1))+1, :);

                plot(geolist{ii}.s, 'FaceColor', faceCol, 'FaceAlpha', 1, 'EdgeColor', 'none');
            end

            % -------- Axes cosmetics --------
            % xlabel(axH, 'Interpreter','none'); ylabel(axV, 'Interpreter','none');
            axis image; xlim(xRange); ylim(yRange);
            colormap(cmap); colorbar; clim([cmin cmax]);
            hold off;
        end

    end

    methods (Static, Access = protected)
        function in = findInGeometry(AXIS1, AXIS2, g)
            arguments
                AXIS1 (:,:) double
                AXIS2 (:,:) double
                g 
            end
            % Assumes geometry objects expose a polyshape in property 's'
            ps = g.s;
            xv = ps.Vertices(:,1);
            yv = ps.Vertices(:,2);
            
            % Logical mask for grid points inside the shape
            in = inpolygon(AXIS1, AXIS2, xv, yv);
        end
    end

    methods (Access = ?Source1D)
        function [eps_cs, mu_cs, axis_cs] = CrossSection(obj, port)
            arguments
                obj
                port Port1D
            end

            ind_axis_port = find(obj.mesh.label==port.dir(2));
            if isempty(ind_axis_port)
                dispError('Device2D:WrongPortDirection', obj.mesh.label(1), obj.mesh.label(2));
            end
            switch ind_axis_port
                case 1
                    [~, ind_start] = min(abs(obj.mesh.axis2.v-port.p1(2)));
                    [~, ind_end] = min(abs(obj.mesh.axis2.v-port.p2(2)));
                    [~, ind_position] = min(abs(obj.mesh.axis1.v-port.p1(1)));

                    eps_cs = obj.eps(ind_position, ind_start:ind_end, :);
                    mu_cs  = obj.mu(ind_position, ind_start:ind_end, :);
                    axis_cs = obj.mesh.axis2.cutAxis(ind_start:ind_end);
                case 2
                    [~, ind_start] = min(abs(obj.mesh.axis1.v-port.p1(1)));
                    [~, ind_end] = min(abs(obj.mesh.axis1.v-port.p2(1)));
                    [~, ind_position] = min(abs(obj.mesh.axis2.v-port.p1(2)));

                    eps_cs = obj.eps(ind_start:ind_end, ind_position, :);
                    mu_cs  = obj.mu(ind_start:ind_end, ind_position, :);
                    axis_cs = obj.mesh.axis1.cutAxis(ind_start:ind_end);
            end
        end
    end
end
