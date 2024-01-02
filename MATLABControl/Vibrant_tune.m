%% tune OPO
%
% comments...

function Vibrant_tune(PA_wavelength)

global Vibrant

if PA_wavelength < 710
    OPOconfig = uint8(0);       % to signal...
else
    OPOconfig = uint8(1);       % ...or idler
end
Vibrant_pointers;
[Vibrant.wllo, Vibrant.wlhi, Vibrant.ierr] = calllib('opotek', 'Config', ...
    OPOconfig, Vibrant.pwllo, Vibrant.pwlhi, Vibrant.pierr);
if (PA_wavelength >= Vibrant.wllo && PA_wavelength <= Vibrant.wlhi)
    % wavelength within range, tune
    action = int8(0);
    par1 = int16(10*PA_wavelength);
    par2 = int16(0);
    Vibrant_pointers;
    Vibrant.ierr = calllib('opotek', 'Motor', action, par1, par2, Vibrant.pierr);
else
    disp( 'Wavelength error; returning')
    return
end
