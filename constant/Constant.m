classdef Constant
    %Constant: SI physical-constant container for solver/material computations
    %
    % Key Properties (SetAccess = private):
    %   v     - Numeric constant value (SI)
    %   name  - Constant identifier: 'c0'|'eta0'|'eps0'|'mu0'
    %   u     - Unit string in internal SI-notation (for display/debug)
    %
    % Key Methods:
    %   Constant(name)  - Construct constant object by identifier and populate v/u

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