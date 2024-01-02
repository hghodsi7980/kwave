
global Vibrant
PAon = 1
PA_wavelength = 680;
PA_Attenuator = 10;
PA_laserenergy = 80; 
PA_RepRate_int = 0;%
%% Opotek Vibrant laser
% #('Vibrant', com#:laser, com#:motors, com#:spectrometer)??
if PAon == 1
    disp('Initializing Vibrant');
    if ~libisloaded('opotek')
        loadlibrary('OPOTEK', @OPOproto, 'alias', 'opotek')
    end
    if PA_wavelength < 600
        disp(' smallest allowed value for wavelength = 600 nm')
        return
    end
    if isempty(Vibrant)
        Vibrant_init;
    end
end

%% set attenuator to specified percentage (0-100)
if PAon == 1 && ~strcmp(PA_Attenuator,'no')
    % PA_Attenuator = 100;
    Attenuator_str = num2str(PA_Attenuator);
    disp(['Setting attenuator to ' Attenuator_str '%']);
    action = int8(5); par1 = int16(PA_Attenuator);par2 = int16(0); Vibrant_pointers;
    [Vibrant.ierr] = calllib('opotek', 'Motor', action, par1, par2, Vibrant.pierr);
end
%% set laser energy
if PAon == 1
    str_laserenergy = num2str(PA_laserenergy);
    disp(['Setting laser energy to ' str_laserenergy '%']);
    action = int8(3); par1 = int16(PA_laserenergy); sts = ''; Vibrant_pointers;
    [sts,Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
end

%% set photoacoustical rep rate
if PAon == 1
    PA_RepRate = 10/(PA_RepRate_int + 1);
    PA_RepRate_str = num2str(PA_RepRate);
    disp(['Setting photoacoustical rep rate to ' PA_RepRate_str ' Hz']);
    action = int8(4); par1 = int16(PA_RepRate_int); sts = ''; Vibrant_pointers;
    [sts,Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
end

%% tune OPO
if PAon == 1
    PA_wavelength_str = num2str(PA_wavelength);
    disp(['Tuning OPO to ' PA_wavelength_str ' nm']);
    Vibrant_tune(PA_wavelength);
end

%% start flashlamp and laser
if PAon == 1
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
end

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

%%
if PAon == 1; Vibrant_tune(PA_wavelengtharray(1,iwl)); end

%%
% stop Vibrant
       
            disp('Stopping Vibrant');
            % stop laser
            action = int8(2); par1 = int16(0); sts = ''; Vibrant_pointers;
            [sts Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
            % stop flashlamp
            action = int8(1); par1 = int16(0); sts = ''; Vibrant_pointers;
            [sts Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
            % % set energy to zero
            %         action = int8(3); par1 = int16(0); sts = ''; Vibrant_pointers;
            %         [sts Vibrant.ierr] = calllib('opotek', 'Laser', action, par1, sts, Vibrant.pierr);
 