classdef Bitmap2D < Layout2D
    methods
        function obj = Bitmap2D(bitmap, unit)
            arguments
                bitmap (:,:)
                unit {mustBeMember(unit, {'m', 'cm', 'mm', 'um', 'nm', 'A'})} = 'm'
            end

            % get vertex coordinates
            boundaries = bwboundaries(bitmap);

            num_coord = numel(boundaries);
            coord_data_stack = cell(num_coord, 1);
            for ii = 1:num_coord
                xy_data_ii = cell2mat(boundaries(ii));
                coord_data_stack(ii) = {xy_data_ii};
            end

            % recreate bitmap in polyshape class
            warning('off', 'MATLAB:polyshape:repairedBySimplify'); % avoid warnings 
            % triggered by hole structures
            xy_data1 = cell2mat(coord_data_stack(1));
            obj@Layout2D(polyshape(xy_data1(:,1), xy_data1(:,2)), unit);

            for jj = 2:num_coord
                xy_data = cell2mat(coord_data_stack(jj));
                polygon_part = Polygon(xy_data(:,1), xy_data(:,2), unit);
                obj = combineGeometry(obj, 'Union', polygon_part);
            end
            
        end
    end

end