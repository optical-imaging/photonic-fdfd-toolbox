classdef (Abstract, Hidden) BC < handle
    %BC: Abstract container class for boundary-condition configuration in solvers
    %
    % Purpose:
    %   Common interface for boundary-condition objects (configuration only).
    %   Numerical enforcement is implemented in solver classes (e.g., FDFD).
    %
    % Abstract Methods (must implement):
    %   setPMLThickness(...)  - Set PML thickness parameter(s)
    %   changeType(...)       - Change boundary type (Perfect/PML/Periodic)

    methods (Abstract)
        setPMLThickness
        changeType
    end
end