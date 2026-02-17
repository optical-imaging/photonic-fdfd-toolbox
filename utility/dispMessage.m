function dispMessage(message_id, varargin)
switch message_id
    case 'Model:ModelSolved'
        disp('Model object has been solved successfully.')

    case 'Model2D:ObjectInitialized'
        disp('2D model has been created.');

    case 'LayerDevice:LayerOverlap'
    overlap_seq = varargin{1};
    msg = "LayerDevice z-overlap detected between device seq:";
    for k = 1:size(overlap_seq,1)
        msg = msg + sprintf(" [%d & %d]", overlap_seq(k,1), overlap_seq(k,2));
    end
    disp(msg);


    case 'Source1D:VariableSetup'
        disp('Source has been set and ready to be used.');
end
end
