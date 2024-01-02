function Vibrant_pointers
% creates pointers to the opotek library functions 
% and gives them the same value as the corresponding
% variable
% the opotek library functions need pointers as input

global Vibrant

if isempty(Vibrant), Vibrant_init; end 

Vibrant.plasercom = libpointer('int32Ptr',Vibrant.lasercom);
Vibrant.pmotorcom = libpointer('int32Ptr',Vibrant.motorcom);
Vibrant.pwllo = libpointer('doublePtr',Vibrant.wllo);
Vibrant.pwlhi = libpointer('doublePtr',Vibrant.wlhi);
Vibrant.pierr = libpointer('int32Ptr',Vibrant.ierr);