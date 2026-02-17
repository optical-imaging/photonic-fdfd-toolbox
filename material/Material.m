classdef Material
    %Material: Linear isotropic material container used by Device meshing
    %
    % Key Properties (SetAccess = private):
    %   eps  - Relative permittivity (scalar)
    %   mu   - Relative permeability (scalar)
    %   lam  - Wavelength associated with database lookup (numeric) or 'NotDefined'
    %
    % Key Methods:
    %   Material(input1,input2)  - Construct from (eps,mu) or (material_name, wavelength)
    %   printInfo()              - Print eps/mu/lam information (debug use)
    %
    % Key Methods (Static, private):
    %   createFromValues(eps,mu)                 - Validate and return eps/mu values
    %   createFromDatabase(material_name,lam)    - Load material data and return eps/mu at wavelength
    %   getMaterialFromDatabase(data,lam)        - Exact-match or interpolate eps/mu from database table

    properties (SetAccess = private)
        eps
        mu
        lam
    end

    methods

        % Constructor
        function obj = Material(input1, input2)
            % Material constructor
            %   Constructs an instance of the Material class using either specified permittivity and permeability values,
            %   or by looking up values from a database based on the material name and wavelength.
            %
            %   Syntax:
            %     obj = Material(input1, input2)
            %
            %   Input:
            %     input1 - Either the relative permittivity value (scalar) or the material name (string).
            %     input2 - Either the relative permeability value (optional, default is 1) or the wavelength (for material name).
            %
            %   Output:
            %     obj - An instance of the Material class.

            arguments
                input1
                input2 = [];
            end
            if isnumeric(input1) && isscalar(input1)
                % First method: Enter eps value, optional mu
                if nargin == 1
                    input2 = 1;
                end
                [eps_val, mu_val] = Material.createFromValues(input1, input2);
                lam_val = 'NotDefined';
            elseif ischar(input1)
                % Second method: Enter material name and wavelength
                [eps_val, mu_val] = Material.createFromDatabase(input1, input2);
                lam_val = input2;
            else
                dispError('Material:WrongInput');
            end

            obj.eps = eps_val;
            obj.mu = mu_val;
            obj.lam = lam_val;
        end

        % Display
        function printInfo(obj)
            info_name = {'Material Name';
                'Epsilon_r';
                'Mu_r';
                'Wavelength'};
            value = {inputname(1);
                num2str(obj.eps);
                num2str(obj.mu);
                [num2str(obj.lam),' m']};

            maxNameLength = max(cellfun(@length, info_name));
            for ii = 1:numel(info_name)
                fprintf('%*s: %-10s\n', maxNameLength+1, info_name{ii}, value{ii});
            end
        end

    end

    methods (Static, Access = private)

        % Method for creating material from eps, optional mu, and label
        function [eps, mu] = createFromValues(eps_val, mu_val)
            arguments
                eps_val (1,1) double
                mu_val (1,1) double = 1
            end

            eps = eps_val;
            mu = mu_val;
        end

        % Method for creating material from material name and wavelength
        function [eps, mu] = createFromDatabase(material_name, wavelength)
            arguments
                material_name char {mustBeMember(material_name, {'Vacuum', 'Si'})}
                wavelength (1,1) double {mustBeReal, mustBePositive}
            end

            % Fetch material properties from the database
            data = load('MaterialData.mat',material_name);
            var_name = fieldnames(data);
            material_data = data.(var_name{1});

            [eps_val, mu_val] = Material.getMaterialFromDatabase(material_data, wavelength);
            eps = eps_val;
            mu = mu_val;
        end

        function [eps, mu] = getMaterialFromDatabase(material_data, wavelength)
            % materialData: Matrix where
            %   - column 1 is wavelength
            %   - column 2 is eps
            %   - column 3 is mu

            wavelengths = material_data(:,1);
            eps_values = material_data(:,2);
            mu_values = material_data(:,3);

            % Check if the wavelength is directly available in the data
            if any(wavelengths == wavelength)
                % Get the index of the exact match
                index = find(wavelengths == wavelength);
                eps = eps_values(index);
                mu = mu_values(index);
            else
                % Check if the wavelength is within the range
                if wavelength >= min(wavelengths) && wavelength <= max(wavelengths)
                    % Use interpolation to estimate eps and mu
                    eps = interp1(wavelengths, eps_values, wavelength);
                    mu = interp1(wavelengths, mu_values, wavelength);
                else
                    dispError('Material:WavelengthOutOfRange');
                end
            end
        end

    end

end