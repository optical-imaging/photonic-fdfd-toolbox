classdef Constant
    % Constant: Physics Constant Class
    %
    % Description:
    %   The `Constant` class provides a convenient way to access commonly used
    %   physical constants, all in SI units. The available constants are:
    %   c0 - Speed of light in vacuum
    %   eta0 - Characteristic impedance of vacuum
    %   eps0 - Permittivity of free space
    %   mu0 - Permeability of free space
    %
    %   Each constant is represented by a value (`v`), its name (`name`), and the unit (`u`).
    %
    % Usage:
    %   constant = Constant('c0');
    %
    % Properties:
    %   v - The numerical value of the constant.
    %   name - The name of the constant.
    %   u - The unit of the constant in SI notation.
    %
    % Methods:
    %   Constant(constant_name) - Constructor method to create an instance of the Constant class.
    %     The `constant_name` must be one of the following: 'c0', 'eta0', 'eps0', 'mu0'.



    properties (SetAccess = private)
        v
        name
        u
    end

    methods

        function obj = Constant(constant_name)
            % Constant constructor
            %   Constructs an instance of the Constant class with a specific constant value.
            %
            %   Syntax:
            %     obj = Constant(constant_name)
            %
            %   Input:
            %     constant_name - A string specifying the desired constant.
            %                    Must be one of: 'c0', 'eta0', 'eps0', 'mu0'
            %
            %   Output:
            %     obj - An instance of the Constant class with properties
            %           corresponding to the selected physical constant.

            arguments
                constant_name {mustBeMember(constant_name, {'c0', 'eta0', 'eps0', 'mu0'})}
            end

            switch constant_name
                case 'c0'
                    obj.v = 2.998e8; 
                    obj.u = '[m]*[s]^[-1]';
                case 'eta0'
                    obj.v = 376.730;
                    obj.u = '[m]^[2]*[kg]*[s]^[-3]*[A]^[-2]'; % ohm
                case 'eps0'
                    obj.v = 8.854e-12;
                    obj.u = '[A]^[2]*[s]^[4]*[kg]^[-1]*[m]^[-3]'; % F/m
                case 'mu0'
                    obj.v = 1.257e-6;
                    obj.u = '[m]*[kg]*[s]^[-2]*[A]^[-2]'; % N/A^2
            end
            obj.name = constant_name;
        end

    end

end