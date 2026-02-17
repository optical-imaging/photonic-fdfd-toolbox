function mustBeEmptyOrBit(s)
if ~isempty(s) && ~strcmp(s, 'bit')
    dispError('Model2D:AddBitmap');
end
end