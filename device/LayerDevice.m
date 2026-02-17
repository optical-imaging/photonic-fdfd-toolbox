classdef LayerDevice < Device2D
    %LayerDevice: Semi-3D layer-stacked device built from 2D layouts extruded along z
    %
    % Key Properties (inherited from Device2D/Device):
    %   seq        - Layer sequence indices (vector when combined); used as layer IDs
    %   geo        - Per-layer 2D geometry (Layout2D objects; cell or array indexed by seq)
    %   mat        - Per-layer material definitions (aligned with seq)
    %   mesh       - Grid3D mesh used for semi-3D voxelization
    %   eps        - 3D relative permittivity map, size Nx×Ny×Nz×3 (Yee components)
    %   mu         - 3D relative permeability map, size Nx×Ny×Nz×3 (Yee components)
    %   meshflag   - True when mesh is assigned (non-empty)
    %
    % Key Properties (LayerDevice-specific):
    %   zmin       - Layer bottom z coordinate(s), scalar or vector aligned with seq
    %   zmax       - Layer top z coordinate(s), scalar or vector aligned with seq
    %
    % Key Properties (Dependent):
    %   t          - Layer thickness(es), t = zmax - zmin
    %   z0         - Layer center(s), z0 = (zmax + zmin)/2
    %
    % Key Methods:
    %   LayerDevice(seq)            - Construct a layer container (inherits Device2D constructor contract)
    %   get.t()                     - Return thickness from zmin/zmax (dependent getter)
    %   get.z0()                    - Return center from zmin/zmax (dependent getter)
    %   setEndPoint('zmin',...,'zmax',...) - Set z-span endpoints for one layer (validates zmin < zmax)
    %   setSpan('Span',...,'z0',...)       - Set z-span by center and thickness
    %   meshDevice(grid3d)          - Build 3D eps/mu by meshing each layer in xy then extruding in z
    %   setBitmap(eps_bit,mu_bit,step) - Override bitmap setter to enforce 3D bitmap + 3-step spacing
    %   dispImg(plot_part,...)      - 3D visualization by extruding polyshapes for each layer
    %
    % Key Methods (Static):
    %   combineDevice(layer1,layer2,...) - Combine layers using Device.combineDevice + concatenate zmin/zmax
    %
    % Key Methods (Private helpers):
    %   getsingleLayer(seq)         - Extract one layer's geo/mat into a standalone LayerDevice
    %   createLayerDevice2D(seq)    - Convert one layer into a Device2D for xy meshing
    %   getGeoBySeq(seq)            - Fetch geometry for a given seq from geo (cell/array)
    %
    % Key Methods (Access = ?Model):
    %   layerOverlap()              - Check neighboring-layer overlaps in z (used by ModelSemi3D assembleDevice)
    %
    % Key Methods (Static, Private):
    %   extrudePolyshape(...)       - Render one polyshape region as a 3D extruded patch (dispImg support)

    properties
        zmin = 0
        zmax = 0
    end

    properties (Dependent)
        t
        z0
    end

    methods
        % constructor
        function obj = LayerDevice(layer_seq)
            arguments
                layer_seq (1,1) {mustBeInteger, mustBePositive}
            end
            obj@Device2D(layer_seq);
        end

        % dependent
        function val = get.t(obj)
            val = obj.zmax - obj.zmin;
        end

        function val = get.z0(obj)
            val = (obj.zmax + obj.zmin)/2;
        end

        % device setup
        function setEndPoint(obj, varargin)
            p = inputParser;
            addParameter(p,'zmin',0,@(x) isnumeric(x) && isscalar(x));
            addParameter(p,'zmax',0,@(x) isnumeric(x) && isscalar(x));
            parse(p,varargin{:});

            if p.Results.zmin >= p.Results.zmax
                error("Invalid z-span: require zmin < zmax.")
            end

            obj.zmin  = p.Results.zmin;
            obj.zmax  = p.Results.zmax;
        end

        function setSpan(obj, varargin)
            p = inputParser;
            addParameter(p,'Span',0,@(x) isnumeric(x) && isscalar(x) && x>0);
            addParameter(p,'z0',0,@(x) isnumeric(x) && isscalar(x));
            parse(p,varargin{:});

            obj.zmin  = p.Results.z0 - p.Results.Span/2;
            obj.zmax  = p.Results.z0 + p.Results.Span/2;
        end

        % device meshing
        function meshDevice(obj, grid, bg_mat)
            arguments
                obj
                grid Grid3D
                bg_mat Material = []
            end

            if isempty(bg_mat)
                eps_bg = 1;
                mu_bg = 1;
            else
                eps_bg = bg_mat.eps;
                mu_bg = bg_mat.mu;
            end

            x = grid.axisx; y = grid.axisy; z = grid.axisz;
            x_ds = x.doublesampleAxis;
            y_ds = y.doublesampleAxis;
            z_ds = z.doublesampleAxis;
            [AXISX, AXISY] = ndgrid(x_ds.v, y_ds.v);

            eps_ds = eps_bg*ones(x_ds.n, y_ds.n, z_ds.n);
            mu_ds  = mu_bg*ones(x_ds.n, y_ds.n, z_ds.n);

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
                xy_in = obj.findInGeometry(AXISX, AXISY, g);

                z_in = false(1,1,z_ds.n);
                zmin_ii = obj.zmin(ii); [~, iz1] = min(abs(z_ds.v - zmin_ii));
                zmax_ii = obj.zmax(ii); [~, iz2] = min(abs(z_ds.v - zmax_ii));
                z_in(1,1,iz1:iz2) = true;

                mask = logical(xy_in.*z_in);

                % assign material value
                eps_val = 1; mu_val  = 1;
                if ii <= numel(obj.mat)
                    m = obj.mat(ii);
                    if isnumeric(m)
                        eps_val = real(m);
                    else
                        if isprop(m,'eps') && ~isempty(m.eps), eps_val = m.eps; end
                        if isprop(m,'mu')  && ~isempty(m.mu),  mu_val  = m.mu;  end
                    end
                end
  
                eps_ds(mask) = eps_val;
                mu_ds(mask)  = mu_val;
            end

            [eps1, eps2, eps3, mu1, mu2, mu3] = obj.sampleMapAll(eps_ds, mu_ds);
            obj.eps = cat(4, eps1, eps2, eps3);
            obj.mu  = cat(4, mu1, mu2, mu3);
            obj.mesh = grid;
        end

        function setBitmap(obj, eps_bit, mu_bit, step_size)
            arguments
                obj
                eps_bit (:,:,:) {mustBeRealAndGE1}
                mu_bit (:,:,:) {mustBeRealAndGE1}
                step_size (1,3) double {mustBeReal, mustBePositive}
            end
            % eps and mu size match
            if ~isequal(size(eps_bit), size(mu_bit))
                dispError('Device:BitmapSizeMismatch');
            end

            obj.eps = repmat(eps_bit, [1 1 1 3]);
            obj.mu  = repmat(mu_bit, [1 1 1 3]);
            obj.mesh = Grid3D(size(eps_bit), step_size);
        end


        % display
        function dispImg(obj, plot_part, varargin)
            arguments
                obj
                plot_part {mustBeMember(plot_part, {'eps','mu'})} = 'eps'
            end
            arguments (Repeating)
                varargin
            end

            p = inputParser;
            addParameter(p,'Colormap',jet(256));
            addParameter(p,'CLim',[1 15],@(x) isnumeric(x) && numel(x)==2);
            addParameter(p,'ShowCaps',true,@(x) islogical(x) && isscalar(x));

            % Subtle edges by default (visible but not dominant)
            addParameter(p,'EdgeColor',[0 0 0]);
            addParameter(p,'LineWidth',0.12,@(x) isnumeric(x) && isscalar(x) && x>=0);
            addParameter(p,'EdgeAlpha',0.10,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);

            % Transparency controls
            addParameter(p,'SideAlpha',0.50,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
            addParameter(p,'CapAlpha',0.90,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1); % key for "occluding" buried parts

            % Visualization controls
            addParameter(p,'View',[35 18],@(v) isnumeric(v) && numel(v)==2); % [az el]
            addParameter(p,'ZScale',6,@(s) isnumeric(s) && isscalar(s) && s>0);
            addParameter(p,'Projection','perspective',@(s) ismember(lower(string(s)), ["orthographic","perspective"]));
            addParameter(p,'SortMethod','depth',@(s) ismember(lower(string(s)), ["depth","childorder"])); % depth helps blending
            addParameter(p,'Aspect','pretty',@(s) ismember(lower(string(s)), ["exact","pretty"]));
            parse(p,varargin{:});

            cmap  = p.Results.Colormap;
            climv = p.Results.CLim;

            figure; hold on; grid on;
            set(gcf,'Renderer','opengl');
            ax = gca;
            set(ax,'SortMethod',char(lower(string(p.Results.SortMethod))));
            set(ax,'Projection',char(lower(string(p.Results.Projection))));

            xlabel('x'); ylabel('y'); zlabel('z');
            switch plot_part
                case 'eps', title('\epsilon_r');
                case 'mu',  title('\mu_r');
            end

            view(p.Results.View(1), p.Results.View(2));
            axis vis3d;

            for ii = 1:numel(obj.seq)
                seq = obj.seq(ii);
                g = obj.getGeoBySeq(seq);
                if isempty(g) || ~isprop(g,'s') || isempty(g.s), continue; end

                % scalar for colormap
                val = 1;
                if ii <= numel(obj.mat)
                    m = obj.mat(ii);
                    if isnumeric(m)
                        val = real(m);
                    else
                        if strcmp(plot_part,'eps') && isprop(m,'eps') && ~isempty(m.eps), val = real(m.eps); end
                        if strcmp(plot_part,'mu')  && isprop(m,'mu')  && ~isempty(m.mu),  val = real(m.mu);  end
                    end
                end

                tt = (val - climv(1)) / (climv(2) - climv(1));
                tt = max(0, min(1, tt));
                faceCol = cmap(round(tt*(size(cmap,1)-1))+1, :);

                z1 = obj.zmin(ii);
                z2 = obj.zmax(ii);

                regs = regions(g.s);
                for rr = 1:numel(regs)
                    LayerDevice.extrudePolyshape(regs(rr), z1, z2, faceCol, ...
                        p.Results.ShowCaps, ...
                        p.Results.EdgeColor, p.Results.LineWidth, p.Results.EdgeAlpha, ...
                        p.Results.SideAlpha, p.Results.CapAlpha);
                end
            end

            colormap(cmap);
            colorbar; clim(climv);

            switch lower(string(p.Results.Aspect))
                case "exact"
                    axis equal
                case "pretty"
                    axis tight
                    daspect([1 1 1/p.Results.ZScale]);
            end

            lighting gouraud;
            camlight headlight;
            camlight(60, 20); % second light helps depth cues
            material dull;

            hold off;
        end

    end

    methods (Static)
        function combined_device = combineDevice(varargin)
            arguments (Repeating)
                varargin LayerDevice
            end
            % thickness info collection
            device_list = [varargin{:}];
            zmin_list  = [device_list.zmin];
            zmax_list  = [device_list.zmax];

            % combine device list and set new one
            combined_device = combineDevice@Device(varargin{:});
            combined_device.zmin = zmin_list;
            combined_device.zmax = zmax_list;
        end
    end

    methods (Access = private)
        function single_layer = getsingleLayer(obj, seq)
            single_layer = LayerDevice(seq);
            if iscell(obj.geo)
                single_layer.geo = obj.geo{seq};
            else
                single_layer.geo = obj.geo(seq);
            end
            single_layer.mat = obj.mat(seq);
        end

        function layerdevice_2d = createLayerDevice2D(obj, seq)
            s = obj.getsingleLayer(seq);
            layerdevice_2d = Device2D(s.seq);
            layerdevice_2d.geo = s.geo;
            layerdevice_2d.mat = s.mat;
        end

        function g = getGeoBySeq(obj, seq)
            if isempty(obj.geo)
                g = [];
            elseif iscell(obj.geo)
                g = obj.geo{seq};
            else
                g = obj.geo(seq);
            end
        end
    end

    methods (Access = ?Model)
        function [isoverlap, overlap_seq] = layerOverlap(obj)
            seq  = obj.seq(:);
            zmin_val = obj.zmin(:);
            zmax_val = obj.zmax(:);

            % enforce ordering within each layer
            zlo = min(zmin_val, zmax_val);
            zhi = max(zmin_val, zmax_val);

            % sort by zmin
            [zlo, idx] = sort(zlo, 'ascend');
            zhi = zhi(idx);
            seq = seq(idx);

            % check overlap only between neighboring layers
            overlap_seq = zeros(0,2);
            for k = 1:numel(seq)-1
                if zhi(k) > zlo(k+1)
                    pair = sort([seq(k), seq(k+1)]);   % <-- added normalization
                    overlap_seq(end+1,:) = pair;       %#ok<AGROW>
                end
            end

            isoverlap = ~isempty(overlap_seq);
        end
    end

    % display
    methods (Static, Access = private)
        function extrudePolyshape(ps, z1, z2, faceCol, showCaps, edgeCol, lw, edgeA, sideA, capA)
            V = ps.Vertices;
            if size(V,1) < 3, return; end

            % Ensure closed loop for side faces
            x = V(:,1); y = V(:,2);
            if x(1) ~= x(end) || y(1) ~= y(end)
                x = [x; x(1)];
                y = [y; y(1)];
            end

            n = numel(x) - 1;
            if n < 3, return; end

            % ---- Side walls ----
            vb = [x(1:end-1), y(1:end-1), z1*ones(n,1)];
            vt = [x(1:end-1), y(1:end-1), z2*ones(n,1)];
            verts_side = [vb; vt];

            i1 = (1:n).';
            i2 = [2:n, 1].';
            faces_side = [i1, i2, i2+n, i1+n];

            % Slightly darken sidewalls to improve depth perception
            sideCol = 0.80 * faceCol;

            patch('Vertices', verts_side, 'Faces', faces_side, ...
                'FaceColor', sideCol, ...
                'EdgeColor', edgeCol, 'LineWidth', lw, ...
                'FaceAlpha', sideA, 'EdgeAlpha', edgeA, ...
                'BackFaceLighting','reverselit', ...
                'FaceLighting','gouraud', ...
                'AmbientStrength',0.25, ...
                'DiffuseStrength',0.80, ...
                'SpecularStrength',0.20, ...
                'SpecularExponent',12);

            if ~showCaps, return; end

            % ---- Caps (top/bottom) ----
            TR = triangulation(ps);
            if isempty(TR) || isempty(TR.ConnectivityList), return; end
            T = TR.ConnectivityList;
            P = TR.Points;

            vb2 = [P, z1*ones(size(P,1),1)];
            vt2 = [P, z2*ones(size(P,1),1)];

            % Make cap edges even more subtle than side edges
            edgeA_cap = 0.6 * edgeA;
            lw_cap    = max(lw*0.8, 0.05);

            % Bottom cap
            patch('Vertices', vb2, 'Faces', T, ...
                'FaceColor', faceCol, ...
                'EdgeColor', edgeCol, 'LineWidth', lw_cap, ...
                'FaceAlpha', capA, 'EdgeAlpha', edgeA_cap, ...
                'BackFaceLighting','reverselit', ...
                'FaceLighting','gouraud', ...
                'AmbientStrength',0.25, ...
                'DiffuseStrength',0.80, ...
                'SpecularStrength',0.20, ...
                'SpecularExponent',12);

            % Top cap (reverse winding so its normal points outward)
            patch('Vertices', vt2, 'Faces', fliplr(T), ...
                'FaceColor', faceCol, ...
                'EdgeColor', edgeCol, 'LineWidth', lw_cap, ...
                'FaceAlpha', capA, 'EdgeAlpha', edgeA_cap, ...
                'BackFaceLighting','reverselit', ...
                'FaceLighting','gouraud', ...
                'AmbientStrength',0.25, ...
                'DiffuseStrength',0.80, ...
                'SpecularStrength',0.20, ...
                'SpecularExponent',12);
        end

    end
end