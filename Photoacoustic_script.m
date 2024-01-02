% This is a modified from verasonic example by Hamed Ghodsi on Dec 2023 for
% using with our photoacoustic setup
 
clear
close all
 
%% === Set commonly modified parameters ========================================
 
P(1).startDepth = 0;          % Acquisition start depth in wavelengths
P(1).endDepth = 128;          % Acquisition end depth may be different for PA
  
% Set PA parameters
flash2Qdelay = 200; % microseconds between trigger input and start of acquisition (which outputs a trigger pulse at time=0)
ne = 1;         % ne = number of acquisitions in PA ensemble for coherent addition in I/Q buffer.
PA_Angle = 0;   % angle of transmit plane wave for use in testing the PA mode in simulation
PA_PRF = 100;   % PA PRF in Hz. To be set in accordance with laser rep rate, when not using the input trigger mode.
                % When using the input trigger mode, remove the TTNA (SeqControl(4)) from the PA events to avoid
                % "missed TTNA" messages.
disp(' *** PhotoAcoustic mode: Using one-way receive-only reconstruction ***') 
% ===============================================================================
 
%% Specify system parameters(Vantage 128)
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 1;    %always start in simulate mode for safety. Hardware mode : 0
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
 
%% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'mm';             % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % computeTrans is used for known transducers.
nElem = Trans.numelements;
Trans.maxHighVoltage = 50;      % set a reasonable high voltage limit.
 
%% Specify PData structure arrays.
% - PA PData structure
PData(1).PDelta(1) = 1.0;
PData(1).PDelta(3) = 0.5;
PData(1).Size(1,1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % PA window rows
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));  % PA window columns
PData(1).Size(1,3) = 1;             % single image page
PData(1).Origin = [-Trans.spacing*63.5,0,P(2).startDepth]; % x,y,z of upper lft crnr.
 
%% Specify Media object and point displacement function
pt1;
Media.function = 'movePoints';
 
%% Specify Resources.
% - RcvBuffer(1) is for PA acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*(na + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for PA reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for PA.
% ImageBuffer(1) is for PA image.
Resource.ImageBuffer(1).datatype = 'double';    % image buffer for PA
Resource.ImageBuffer(1).numFrames = Resource.ImageBuffer(1).numFrames;
% DisplayWindow
Resource.DisplayWindow(1).Title = mfilename;
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;
 
% ------Specify structures used in Events------
%% Specify Transmit waveforms structure
% - PA transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,6,1];
 
%% Specify Transmit beams structure
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 1); % na TXs for 2D + 1 for PA
% -- only one TX struct needed for PA
TX(1).waveform = 1;
TX(1).Apod = zeros(1,Trans.numelements);     % THIS COMMAND TURNS OFF ALL TRANSMITTERS because we are in PA mode
 
%% Specify TPC structures.
% This allows one to use different transmit profile for PA ... only relevant if transmitters are active
% --- TPC(1) is not used ---
TPC(1).name = 'PA';
TPC(1).maxHighVoltage = 35;
 
%% Analog front end gain settings.
RcvProfile(1).LnaGain = 24;     % 12, 18, or 24 dB 
RcvProfile(1).condition = 'immediate';
 
%% Specify Receive structure arrays.
maxAcqLngthPA =  sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(1).startDepth;
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128), ...
                        'TGC', 1, ...           % TGC(1) is tied to the GUI sliders
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, Resource.RcvBuffer(1).numFrames);
 
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = ne*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1; % move points before doing ensemble of different angle plane waves
    % PA acquisitions
    for j = 1:ne
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % There could be multiple acqusitions per frame
        Receive(j+k).startDepth = P(1).startDepth;
        Receive(j+k).endDepth = P(1).startDepth + wl4sPer128*ceil(maxAcqLngthPA/wl4sPer128);
        Receive(j+k).TGC = 1;         
    end
end
 
%% Specify TGC Waveform structures.
% - PA TGC
TGC(1).CntrlPts =[ 0 272 662 662 662 662 662 662]; % TO BE MODIFIED HERE AFTER COLLECTING PA DATA
TGC(1).rangeMax = P(2).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(2));
 
%% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for PA. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', zeros(1,1)), 1, 2);
% - Set Recon values for PA ensemble.
Recon(1).pdatanum = 1;
Recon(1).IntBufDest = [1,1];
Recon(1).ImgBufDest = [1,-1];
Recon(1).RINums(1,1:ne) = 1:ne;   % 'ne' ReconInfos needed for PA ensemble.
 
%% Define ReconInfo structures.
% - For PA, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % 4=accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, ne);
% - Set specific ReconInfo attributes.
%  - ReconInfos for PA ensemble.
for j = 1:ne
    if j==1, ReconInfo(j).mode = 'replaceIQ'; end
    ReconInfo(j).txnum = na + 1;
    ReconInfo(j).rcvnum = na + j;
end
if ne>1
    ReconInfo(ne).mode = 'accumIQ_replaceIntensity'; % 5=accum and detect
else
    ReconInfo(ne).mode = 'replaceIntensity'; % 0=replace IQ data, detect, and replace Intensity data;  1=Add the new reconstructed intensity data to the data in the ImageBuffer
end
 
%% Specify Process structure arrays.
cpt = 22;       % define here so we can use in UIControl below
cpers = 80;     % define here so we can use in UIControl below
 
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',10,...      % reject level
                         'persistMethod','dynamic',...
                         'persistLevel',cpers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','upperHalf',...
                         'threshold',cpt,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
 
 
%% Specify SeqControl structure arrays.
% -- Change to Profile (PA)
SeqControl(1).command = 'setTPCProfile';
SeqControl(1).condition = 'next';
SeqControl(1).argument = 1;
% output trigger
SeqControl(2).command = 'triggerOut';
% -- Jump back to start.
SeqControl(3).command = 'jump';
SeqControl(3).argument = 1;
% set receive profile
SeqControl(4).command = 'setRcvProfile';
SeqControl(4).argument = 1;
% -- PRF for PA ensemble
SeqControl(5).command = 'timeToNextAcq';
SeqControl(5).argument = round(1/(PA_PRF*1e-06)); % (10 msecs for PA_PRF=100 Hz)
% input trigger
SeqControl(6).command = 'triggerIn';
SeqControl(6).condition = 'Trigger_1_Rising'; % Trigger input 1, enable with rising edge
SeqControl(6).argument = 2; % 500 msec timeout delay
    % (Timeout range is 1:255 in 250 msec steps; 0 means timeout disabled)
% noop delay between trigger in and start of acquisition
SeqControl(7).command = 'noop';
SeqControl(7).argument = fix(flash2Qdelay)*5; % noop counts are in 0.2 microsec increments
% sync command
SeqControl(8).command = 'sync';
 
% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.
nsc = 9;  % next SeqControl number
 
% Specify factor for converting sequenceRate to frameRate. Every 3 frames
% exit to matlab
frameRateFactor = 3;
 
%% Specify Event structure arrays.
n = 1;  % event number 
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire PA ensemble.
    for j = 1:ne
        % set a triger out event to the laser
        Event(n).info = 'Trigger OUT event';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        % Wait for input trigger from flash lamp firing
        Event(n).info = 'Wait for Trigger IN';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 6;
        n = n+1;
 
        % Pause for optical buildup
        Event(n).info = 'noop and sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [7,8];
        n = n+1;
 
        % send trigger output at start of every PA acquisition to fire Q-switch
        Event(n).info = 'Acquire PA event';
        Event(n).tx = 1;
        Event(n).rcv = ne*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 5;
        n = n+1;
    end
 
    Event(n).info = 'Transfer Data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    n = n+1;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer (each needs a different value of nsc)
      nsc = nsc+1;
 
    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;
 
    Event(n).info = 'PA image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab only every 3rd frame to prevent slowdown
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end
 
Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;
 
 
%% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%-UI#1Callback');
 
% - Color Priority Threshold Slider
UI(2).Control = {'UserB2','Style','VsSlider','Label','Color Priority','SliderMinMaxVal',[0,255,cpt],...
                 'SliderStep',[1/255,0.1],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%-UI#2Callback');
 
% - Color Persistence Slider
UI(3).Control = {'UserB1','Style','VsSlider','Label','Color Persistence','SliderMinMaxVal',[0,100,cpers],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%-UI#3Callback');
 
 
%% Save all the structures to a .mat file, and run VSX automatically
filename = ('L11-5vFlashPA');   % define variable 'filename' to permit VSX to skip user query for matfile
save (['MatFiles/',filename])                 % save the structures to a matfile
% VSX                             % invoke VSX automatically when running this Setup script
 
return
 
 
%% **** Callback routines to be encoded by text2cell function. ****
% ---------------------------------------------------------------------
%-UI#1Callback - Sensitivity cutoff change
    ReconL = evalin('base', 'Recon');
    for i = 1:size(ReconL,2)
        ReconL(i).senscutoff = UIValue;
    end
    assignin('base','Recon',ReconL);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon'};
    assignin('base','Control', Control);
    return
%-UI#1Callback
 
 
%-UI#2Callback - Color Threshold change
    % Set the value in the Process structure for use in cineloop playback.
    Process = evalin('base','Process');
    for k = 1:2:length(Process(1).Parameters)
        if strcmp(Process(1).Parameters{k},'threshold'), Process(1).Parameters{k+1} = UIValue; end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.threshold.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',2,'threshold',UIValue};
    assignin('base','Control', Control);
%-UI#2Callback
 
%-UI#3Callback - Color Persistence change
    % Set the value in the Process structure for use in cineloop playback.
    Process = evalin('base','Process');
    for k = 1:2:length(Process(1).Parameters)
        if strcmp(Process(1).Parameters{k},'persistLevel'), Process(1).Parameters{k+1} = UIValue; end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',1,'persistLevel',UIValue};
    assignin('base','Control', Control);
%-UI#3Callback

