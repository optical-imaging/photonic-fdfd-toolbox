classdef GDSII < Layout2D
    %GDSII: Planar geometry imported from GDSII file
    %
    % Purpose:
    %   Wraps GDSII polygon data as a single planar geometry.
    %
    % Key Properties:
    %   s  - polyshape representing imported geometry
    %
    % Key Methods:
    %   GDSII(filename)         - Load GDSII file and build geometry

    methods
        function obj = GDSII(GDSII_file_name)
            arguments
                GDSII_file_name char
            end

            % load .gds
            library = read_gds_library(GDSII_file_name); % load .gds, get gds_library
            structure = get(library, 1);  % get gds_structure from gds_library
            element = get(structure);

            % extract polygon coordinates
            num_element = numel(element);
            coord_data_stack = cell(num_element, 1);
            length_unit = get(library, 'uunit'); % convert length unit into [m]
            for ii = 1:num_element
                xy_data_ii = cell2mat(get(element{ii}, 'xy')).*length_unit;
                coord_data_stack(ii) = {xy_data_ii};
            end

            % recreate polygon in MATLAB
            warning('off', 'MATLAB:polyshape:repairedBySimplify'); % avoid warnings
            % triggered by hole structures
            xy_data1 = cell2mat(coord_data_stack(1));
            obj@Layout2D(polyshape(xy_data1(:,1), xy_data1(:,2)));

            for jj = 2:num_element
                xy_data = cell2mat(coord_data_stack(jj));
                polygon_part = Polygon(xy_data(:,1), xy_data(:,2));
                obj = combineGeometry(obj, 'Union', polygon_part);
            end
        end
    end

end