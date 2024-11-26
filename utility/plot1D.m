function plot1D(axis_plot, curve, varargin)
arguments
    axis_plot Axis
    curve (:,1) double
end

arguments (Repeating)
    varargin
end

p = inputParser;
addParameter(p, 'DataType', 'real', ...
    @(x) ismember(x,{'real','imag','abs'})); % Default plot is real part
addParameter(p, 'Unit', axis_plot.u, ...
    @(x) ismember(x,{'m', 'cm', 'mm', 'um', 'nm', 'A'})) % Unit of the first axis is set as default
addParameter(p, 'Colorbar', 'off', ...
    @(x) ismember(x,{'on','off'})); % turn off colorbar as default
parse(p, varargin{:});

data_mode = p.Results.DataType;
unit_info = p.Results.Unit;
colorbar_flag = p.Results.Colorbar;

axis_plot = axis_plot.changeUnit(unit_info);


switch data_mode
    case 'real'
        plot(axis_plot.v, real(curve), '-', 'LineWidth', 1.5);
    case 'imag'
        plot(axis_plot.v, imag(curve), '-', 'LineWidth', 1.5);
    case 'abs'
        plot(axis_plot.v, abs(curve), '-', 'LineWidth', 1.5);
end

xlabel([axis_plot.label, ' (Unit:', unit_info, ')']);
axis tight;
if strcmp(colorbar_flag, 'on')
    colorbar;
end
end