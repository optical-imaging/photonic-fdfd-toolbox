function plot2D(mesh, image, varargin)
arguments
    mesh Mesh2D
    image (:,:) double
end

arguments (Repeating)
    varargin
end

p = inputParser;
addParameter(p, 'DataType', 'real', ...
    @(x) ismember(x,{'real','imag','abs'})); % Default plot is real part
addParameter(p, 'Colorbar', 'off', @(x) ismember(x,{'on','off'}))% Default colorbar is off
addParameter(p, 'Unit', mesh.axis1.u, ...
    @(x) ismember(x,{'m', 'cm', 'mm', 'um', 'nm', 'A'})) % Unit of the first axis is set as default
parse(p, varargin{:});

data_mode = p.Results.DataType;
colorbar_flag = p.Results.Colorbar;
unit_info = p.Results.Unit;

axis1 = mesh.axis1.changeUnit(unit_info);
axis2 = mesh.axis2.changeUnit(unit_info);

switch data_mode
    case 'real'
        imagesc(axis1.v, axis2.v, real(image)');
    case 'imag'
        imagesc(axis1.v, axis2.v, imag(image)');
    case 'abs'
        imagesc(axis1.v, axis2.v, abs(image)');
end

xlabel([mesh.axis1.label, ' (Unit:', unit_info, ')']);
ylabel(mesh.axis2.label);
axis image;
colormap jet;
if strcmp(colorbar_flag, 'on')
    colorbar;
end
end