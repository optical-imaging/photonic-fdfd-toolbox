function dispError(error_id, varargin)
switch error_id
    case 'Model2D:ChangeLabel'
        error(error_id, 'Plane Label can not be changed after device has been assigned.');
    case 'Model2D:ChangeWavelength'
        error(error_id, 'Wavelength can not be changed after device has been assigned.');
    case 'Model2D:AddDevice'
        error(error_id, 'Before assigning device, plane label and wavelength must be set first.');
    case 'Model2D:SetMesh'
        error(error_id, 'Before assigning mesh, plane label, wavelength, and device must be set first.');
    case 'Model2D:SetSolver'
        error(error_id, 'Before assigning solver, plane label, wavelength, device, and mesh must be set first.');
    case 'Model2D:AddPort'
        error(error_id, 'Before assigning port, plane label, wavelength, device, mesh, and solver must be set first.');
    case 'Model2D:NoPortNeed'
        error(error_id, 'There is no port in EigenMode simulation.');
    case 'Model2D:NoSourceNeed'
        error(error_id, 'There is no source in EigenMode simulation.');
    case 'Model2D:Solve'
        error(error_id, 'Please set all parameters first before run the solver.');
    case 'Model2D:ScatteringNotSolved'
        error(error_id, 'Scattering problem has not been solved. Please solve model first.');
    case 'Model2D:EigenModeNotSolved'
        error(error_id, 'Eigenmode problem has not been solved. Please solve model first.');

    case 'Geometry1D:EndpointsOverlap'
        error(error_id, 'The coordinates of the two ends of the line overlap.');
    case 'Geometry1D:NotParallel'
        error(error_id, ['The No.', varargin{1}, ' segment line is not parallel to the first one.']);
    case 'Geometry1D:NotCollinear'
        error(error_id, ['The No.', varargin{1}, ' segment line is not collinear to the first one.']);
    case 'Port1D:RangeSidesEqual'
        error(error_id, 'Maximum and minimum values of port range have to be different.');
    case 'Port1D:DirectionWrong'
        error(error_id, ['Direction of port should be within the Plane. Please choose ', '''', varargin{1}, '''', ' or ', '''', varargin{2}, '''', '.']);
    case 'Ring:RingRadiusEqual'
        error(error_id, 'Radius of inner and outer circles have to be different.');

    case 'Material:WrongInput'
        error(error_id,['Invalid input. Expected numeric eps or mu value,' ...
            ' or material name and wavelength.']);
    case 'Material:WavelengthOutOfRange'
        error(error_id, ['The wavelength is out of validated wavelength' ...
            ' range of this material']);
        
    case 'Device2D:MeshDevice'
        error(error_id, 'Please mesh device first before going into FDFD solver.');

    case 'Mesh2D:SameAxisLabel'
        error(error_id,'Please input axises along two different dimensions.');
    case 'Mesh2D:AxisNotIncluded'
        error(error_id,['This 2D region is difined on ', varargin{1},...
            ' plane. Please select either ', varargin{2}, ' or ',varargin{3},' axis.']);
    case 'Mesh3D:AxisMistake'
        error(['The three axises constructing a Mesh3D object have to be' ...
            ' x, y, and z axis.']);
    
    case 'Device:LayerOverlap'
        error(error_id, 'Devices to be combined have to be on different layers.');
    case 'Device1D:MeshNotMatch'
        error(error_id,['Device ','''',varargin{1},'''',' has to be defined' ...
            ' on the same mesh as previous device.']);

    case 'FDFD2D:DefineSourceFirst'
        error(error_id, 'Please define the source object before placing it in FDFD solver.');
end
end