function dispWarning(warning_id, varargin)
switch warning_id
    case 'Geometry2D:UnitsNotAlign'
        warning(warning_id, ['Units of input geometries are different. They are unified to ','''',varargin{1},'''','.']);

    case 'Mesh2D:DifferentAxisUnit'
        warning(warning_id, ['Since two axes of the mesh have different unit settings' ...
            'the area is calculated with the unit of the first axis as default.']);
    case 'Mesh3D:DifferentAxisUnit'
        warning(warning_id, ['Since three axes of the mesh have different unit settings' ...
            'the volume is calculated with the unit of the first axis as default.']);
end
end