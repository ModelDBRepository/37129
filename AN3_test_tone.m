% Simple Test for AN3 code     10/19/03  LHC
% Response of 1 kHz CF fiber to a 100 Hz tone
clear all
close all

CF = 1000; % Hz
tonefreq = 1000;  % Hz

srate = 50000; % sampling rate (Hz)
tdres = 1/srate;
dur = 0.050; % stimulus duration in seconds 
ramplength = 0.010; % ramp duration in seconds (note - duration includes 1/2 of ramp)
rampts = ramplength * srate;
steadypts = dur * srate - rampts;
totalpts = steadypts + (rampts*2);
% WAVEFORM ENVELOPE - Cosine^2 On/off ramps
step = pi/(rampts-1);
x=[0:step:pi];
offramp = (1+cos(x))./2;
onramp = (1+cos(fliplr(x)))./2;
o=ones(1,steadypts);
wholeramp = [onramp o offramp]; % Envelope for stimulus (i.e. on/off ramps)
quiet = zeros(1,floor(0.050/tdres)); % use this to add silence at end of stimulus
clear step x;

plotnum = 0;
for spl = 0:20:80 
    plotnum = plotnum + 1;
    % Amplitude scaling of tone (into Pascals for AN model)
    aa = (1/0.707) * 20.e-6 * 10.0^(spl / 20.0); % Scale tone into Pascals <- Key !!!
    step = tonefreq * 2.0 * pi * tdres;
    x=[0:step:((totalpts-1)*step)]; % This is a speedy way to create a tone.
    inputsound = aa * sin(x); %scale tone 
    inputsound = [(inputsound .* wholeramp) quiet] ; % mult by on/off ramps, and add quiet period at end
    clear step x;
    
    sout=anmod3m(CF, inputsound);

    subplot(5,1,plotnum)
    plot((1:length(sout))*tdres*1e3,sout)
    set(gca,'Fontsize',14)
    xlabel('Time (msec)')
    ylabel('Sp/sec')
    text(1.25 * dur * 1e3,max(sout)/2.,[num2str(spl) 'dB SPL']);
    if plotnum == 1
        title(['CF=' num2str(CF) ' Hz - Response to Tone at CF'])
    end

end