function Vibrant_init
% loads opotek library and initializes the Vibrant
global Vibrant

if ~libisloaded('opotek')
    loadlibrary('OPOTEK', @OPOproto64, 'alias', 'opotek')
end

wllo = -99; wlhi = -99; %
motorcom = int32(0); 
lasercom = int32(0);
ierr = int32(0);
Vibrant.inifile = 'opotek_27043_S02_SN2014_1118_20130227.ini';
Vibrant.motorcom = motorcom;
Vibrant.lasercom = lasercom;
Vibrant.wllo = wllo;
Vibrant.wlhi = wlhi;
Vibrant.ierr = ierr;

Vibrant_pointers; % creates and gives values to pointers

[Vibrant.inifile, Vibrant.motorcom, Vibrant.lasercom, Vibrant.ierr] = ...
    calllib('opotek', 'Init', Vibrant.inifile, ...
    Vibrant.pmotorcom, Vibrant.plasercom, Vibrant.pierr);
