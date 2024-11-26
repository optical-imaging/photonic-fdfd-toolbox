function [unit_factor, unit_name] = LengthUnit(unit_name)
arguments
    unit_name {mustBeMember(unit_name, {'m', 'cm', 'mm', 'um', 'nm', 'A'})}
end
name_list = {'m', 'cm', 'mm', 'um', 'nm', 'A'};
factor_list = [1, 1e-2, 1e-3, 1e-6, 1e-9, 1e-10];
unit_factor = factor_list(strcmp(name_list, unit_name));
end