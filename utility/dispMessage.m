function dispMessage(message_id,varargin)
switch message_id
    case 'Model2D:ObjectInitialized'
        disp('2D model has been created.');
    case 'Source1D:VariableSetup'
        disp('Source has been set and ready to be used.');
end
end