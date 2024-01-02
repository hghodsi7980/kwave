addpath('c:\opotek');
addpath('C:\Program Files (x86)\OPOTEK_Control_Software');
%addpath('D:\LaserControl'); 

global Vibrant
PAon= 1;
PA_wavelength = 650;
%%
%% Opotek Vibrant laser
% #('Vibrant', com#:laser, com#:motors, com#:spectrometer)??
if PAon == 1
    disp('Initializing Vibrant');
    if ~libisloaded('opotek')
        loadlibrary('OPOTEK', @OPOproto64, 'alias', 'opotek')
    end
    if PA_wavelength < 600
        disp(' smallest allowed value for wavelength = 600 nm')
        return
    end
    if isempty(Vibrant)
        Vibrant_init;
    end
end

%% set laser energy
PA_laserenergy = 100;
str_laserenergy = num2str(PA_laserenergy);
disp(['Setting laser energy to ' str_laserenergy '%']);
action = int8(3); par1 = int16(PA_laserenergy); sts = ''; Vibrant_pointers;
[sts,Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);

%% tune OPO
PA_wavelength = 650;
PA_wavelength_str = num2str(PA_wavelength);
disp(['Tuning OPO to ' PA_wavelength_str ' nm']);
Vibrant_tune(PA_wavelength);

%% start flashlamp and laser
%if PAon == 1
    % start flashlamp
    disp('Starting flashlamp');
    disp('ATTENTION: the laser will be turned on in 8 seconds');
    action = int8(1); par1 = int16(1); sts = ''; Vibrant_pointers;
    [sts,Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
    t1 = cputime;
    % incorporate fl.ash lamp stabilization delay
    while cputime - t1 < 8.1  % in seconds
        % wait
    end
    % start laser
    action = int8(2); par1 = int16(1); sts = ''; Vibrant_pointers;
    [sts,Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
%end


%% stop Vibrant
if PAon==1
    disp('Stopping Vibrant');
    % stop laser
    action = int8(2); par1 = int16(0); sts = ''; Vibrant_pointers;
    [sts Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
    % stop flashlamp
    action = int8(1); par1 = int16(0); sts = ''; Vibrant_pointers;
    [sts Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
    % % set energy to zero
    % action = int8(3); par1 = int16(0); sts = ''; Vibrant_pointers;
    % [sts. Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
end