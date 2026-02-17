classdef (Abstract, Hidden) Layout2D < Layout
    %Layout2D: Base class for 2D planar geometry represented by polyshape
    %
    % Purpose:
    %   Wraps MATLAB polyshape for planar device geometry.
    %   Provides boolean geometry ops and cross-section extraction.
    %
    % Key Properties (protected set):
    %   s    - polyshape representing the 2D shape
    %
    % Key Methods:
    %   Layout2D(poly_shape)                 - Construct from a polyshape
    %   combineGeometry(operation, varargin) - Boolean ops: Union/Intersect/Subtract/XOR

    properties (SetAccess = protected)
        s
    end

    methods
        % constructor
        function obj = Layout2D(poly_shape)
            arguments
                poly_shape polyshape
            end

            obj.s = poly_shape;
        end

        % combine polyshapes
        function obj = combineGeometry(obj, operation, varargin)
            arguments
                obj
                operation {mustBeMember(operation, {'Union', 'Intersect', 'Subtract', 'ExclusiveOR'})}
            end
            arguments (Repeating)
                varargin
            end

            switch operation
                case 'Union'
                    for ii = 1:numel(varargin)
                        obj.s = union(obj.s, varargin{ii}.s);
                    end
                case 'Intersect'
                    for ii = 1:numel(varargin)
                        obj.s = intersect(obj.s, varargin{ii}.s);
                    end
                case 'Subtract'
                    for ii = 1:numel(varargin)
                        obj.s = subtract(obj.s, varargin{ii}.s);
                    end
                case 'ExclusiveOR'
                    for ii = 1:numel(varargin)
                        obj.s = xor(obj.s, varargin{ii}.s);
                    end
            end
        end
    end

end