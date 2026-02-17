function dispError(error_id, varargin)
switch error_id
    case 'Model:SetWavelength'
        error(error_id, 'Set wavelength before configuring solver.')
    case 'Model:ChangeWavelength'
        error(error_id, 'Wavelength can not be changed after solver has been set.');
    case 'Model:ChangeBackground'
        error(error_id, 'Background can not be changed after device has been assigned.');
    case 'Model:AddDevice'
        error(error_id, 'Wavelength must be set first before adding device.');
    case 'Model:DeviceMissing'
        error(error_id, 'Selected Device is not existing. Please create it first or check the Device index.');
    case 'Model:NoDevice'
        error(error_id, 'There is no Device for assembling. Please addDevice first.');
    case 'Model:DeviceIndexGap'
        error(error_id, 'Please addDevice in order (1,2,3,...).');
    case 'Model:SetMesh'
        error(error_id, 'Wavelength and device must be set first before setting mesh grid.');
    case 'Model:NoMesh'
        error(error_id, 'No mesh has been set. Please define mesh grid first.');
    case 'Model:MeshDevice'
        error(error_id, 'Mesh device only after setting or editting mesh grid.')
    case 'Model:EigenModeNoPML'
        error(error_id, 'Can''t assign PML boundary to EigenMode solver.');
    case 'Model:AddPort'
        error(error_id, 'Wavelength, device, mesh grid, and solver must be set first before assigning port.');

    case 'Model2D:ChangeLabel'
        error(error_id, 'Plane Label can not be changed after device has been assigned.');
    case 'Model2D:SetSolver'
        error(error_id, 'Plane label, wavelength, device, and mesh must be set first before choosing solvers.');
    case 'Model2D:SetScatteringPol'
        error(error_id, 'For Scattering solver, the 3rd argument must be polarization (''TE'' or ''TM'').');
    case 'Model2D:SetBC'
        error(error_id, 'Before setting boundary condition, plane label, wavelength, device, mesh, and solver must be set first.');
    case 'Model2D:BCaxisWrong'
        error(error_id, ['Please choose ', '''', varargin{1}, '''', ' or ', '''', varargin{2}, '''', ' to set boundary condition.']);
    case 'Model2D:SetAllBC'
        error(error_id, 'Please set boundary condition for all boundaries.');
    case 'Model2D:NoPortNeed'
        error(error_id, 'There is no port in EigenMode simulation.');
    case 'Model2D:NoSourceNeed'
        error(error_id, 'There is no source in EigenMode simulation.');
    case 'Model2D:AddSource'
        error(error_id, 'Add source only after port has been set up.');
    case 'Model2D:SetBeforeSolving'
        error(error_id, 'Please set all parameters first before run the solver.');
    case 'Model2D:ScatteringNotSolved'
        error(error_id, 'Scattering problem has not been solved. Please solve model first.');
    case 'Model2D:EigenModeNotSolved'
        error(error_id, 'Eigenmode problem has not been solved. Please solve model first.');
    case 'Model2D:NoSolutionData'
        error(error_id, 'No solution data. Please solve model first.');

    case 'Modelvar3D:SetSolver'
        error(error_id, 'Plane label, wavelength, device, and mesh must be set first before choosing solvers.');
    case 'Modelvar3D:SetBC'
        error(error_id, 'Before setting boundary condition, plane label, wavelength, device, mesh, and solver must be set first.');
    case 'Modelvar3D:AddSource'
        error(error_id, 'Add source only after port has been set up.');
    case 'Modelvar3D:SetBeforeSolving'
        error(error_id, 'Please set all parameters first before run the solver.');
    case 'Modelvar3D:NoSolutionData'
        error(error_id, 'No solution data. Please solve model first.');

    case 'Scattering2D:BC2DLabelNotMatch'
        error(error_id, 'Labels of input BC2D object and current Scattering2D object don''t match.');
    case 'varFDFD:Mesh2DNotXY'
        error(error_id, 'Labels of input Grid2D object has to be defined on x-y plane for varFDFD solver.');
    case 'varFDFD:BC2DNotXY'
        error(error_id, 'Labels of input BC2D object has to be defined on x-y plane for varFDFD solver.');
    case 'EigenMode1D:BCLabelNotMatch'
        error(error_id, 'Labels of input BC1D object and current EigenMode1D object don''t match.');
    case 'EigenMode2D:BCLabelNotMatch'
        error(error_id, 'Labels of input BC2D object and current EigenMode2D object don''t match.');

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
    case 'Port1D:PortLableNotMatch'
        error(error_id, ['Direction dimension (',varargin{1},') is not included in working plane label (',varargin{2},').']);
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
    case 'Device2D:InputBit'
        error(error_id, "Input bit map must be a logical array.")
    case 'Grid2D:ReorderedLabel'
        error(erro_id, 'Invalid input plane label. Select ''xy'', ''yz'', or ''zx''.');S
    case 'Grid2D:SameAxisLabel'
        error(error_id,'Please input axises along two different dimensions.');
    case 'Grid2D:AxisNotIncluded'
        error(error_id,['This 2D region is difined on ', varargin{1},...
            ' plane. Please select either ', varargin{2}, ' or ',varargin{3},' axis.']);
    case 'Grid3D:AxisMistake'
        error(error_id,['The three axises constructing a Mesh3D object have to be' ...
            ' x, y, and z axis.']);

    case 'BC1D:PMLSetOnly'
        error(error_id, 'Number of layers of PML can only be set when boudanry''s type is ''PML''.');
    case 'BC2D:DuplicateLabel'
        error(error_id, 'The boudary condition settings are duplicated.');
    case 'BC2D:InvalidLabel'
        error(error_id, 'Input label is not included in current BC2D object.')

    case 'Source:PortLabelNotMatch'
        error(error_id, "Working-plane label of port doesn't match with source's working-plane.");

    case 'Device:SeqOverlap'
        error(error_id, 'Devices to be combined have to be with different sequence numbers.');
    case 'Device:BitmapNotEqual'
        error(error_id, 'Sizes of eps and mu bitmaps are not equal.');
    case 'Device:BitmapDimTooMany'
        error(error_id, 'Input bitmap array has to be 1-D, 2-D, or 3-D double array.')
    case 'Device1D:MeshNotMatch'
        error(error_id,['Device ','''',varargin{1},'''',' has to be defined' ...
            ' on the same mesh as previous device.']);
    case 'Device1D:BitmapColumn'
        error(error_id, 'The input bitmap can be either Nx3 (defined as Yee grid) or Nx1 (it will be automatically repmatted as 3 columns).')

    case 'Scattering:SetSourceFirst'
        error(error_id, 'Complete setting Source object before assigning to Scattering object.')
    case 'Scattering2D:ScatteringNotSetUp'
        error(error_id, 'Scattering2D solve hasn''t been set up completely.')
    case 'varFDFD:setvarFDFDFirst'
        error(error_id, 'Please set reference profile and polarization of varFDFD solve first.')
end
end