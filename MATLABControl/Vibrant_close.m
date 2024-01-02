function Vibrant_close

global Vibrant

if ~libisloaded('opotek')
    loadlibrary('OPOTEK', @OPOproto, 'alias', 'opotek')
end
[Vibrant.pierr] = calllib('opotek', 'Close', Vibrant.pierr);
unloadlibrary('opotek')
