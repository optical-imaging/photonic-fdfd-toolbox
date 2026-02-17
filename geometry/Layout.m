classdef (Abstract, Hidden) Layout
    %Layout: Abstract base class for geometry containers used by Device objects
    %
    % Purpose:
    %   Defines the minimal geometry interface shared by all layout types.
    %   Layout objects represent geometry only (no materials, no mesh).
    %
    % Key Methods (Abstract):
    %   combineGeometry(...)  - Combine geometry objects (implementation depends on dimension)

    methods (Abstract)
        combineGeometry
    end

end