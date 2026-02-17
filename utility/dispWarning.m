function dispWarning(warning_id, varargin)
switch warning_id
    case 'Model:PerfectForScattering'
        warning(warning_id, '"Perfect" boudanry condition may cause error in Scattering problem.');
    case 'Model:WavelengthChange'
        warning(warning_id, 'Wavelength has been changed. Please reset solver, port and source, and resolve model.');

    case 'Geometry2D:UnitsNotAlign'
        warning(warning_id, ['Units of input geometries are different. They are unified to ','''',varargin{1},'''','.']);

    case 'Grid2D:DefaultLabel'
        warning(warning_id, 'No grid label specified. Defaulting to ''xy'' (axis1=x, axis2=y).');
    case 'Grid2D:ReorderedLabel'
        warning(warning_id, 'Input label "%s" reordered to "%s" (axis1=%s, axis2=%s).', varargin{:});
    case 'Grid2D:DifferentAxisUnit'
        warning(warning_id, ['Since two axes of the mesh have different unit settings' ...
            'the area is calculated with the unit of the first axis as default.']);
    case 'Grid3D:DifferentAxisUnit'
        warning(warning_id, ['Since three axes of the mesh have different unit settings' ...
            'the volume is calculated with the unit of the first axis as default.']);

    case 'Scattering:SourceLamNotMatch'
        warning(warning_id, ['Wavelength of input Source object is different from wavelength of current Scattering object set before....' ...
            ' Wavelength of Scattering object will be replaced as that of input Source object.'])
end
end