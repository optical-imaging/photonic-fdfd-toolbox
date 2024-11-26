classdef (Abstract) Grid
     % Grid: Even-distributed mesh grid
    %
    % Description:
    %   The `Grid` class is an abstract base class representing an evenly distributed mesh grid.
    %   This class provides an interface for operations related to mesh grids and must be subclassed.
    %
    % Abstract Methods:
    %   changeUnit - Abstract method to change the length unit of the grid.
    %   printInfo - Abstract method to print information about the grid.

    methods (Abstract)
        changeUnit
        printInfo
    end

end