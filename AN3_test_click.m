% Simple Test for AN3 code     10/19/03  LHC
% Response of 1 kHz CF fiber to 100 usec clicks
clear all
close all

CF = 1000; % Hz

srate = 50000; % sampling rate (Hz)
tdres = 1/srate;
dur = 0.025; % stimulus duration in seconds 

iplot = 0;
for spl = 40:20:120 
    iplot = iplot + 1;
    % Amplitude scaling of tone (into Pascals for AN model)
    aa = 20e-6 * 10.0^(spl / 20.0); % Scale peak of click into Pascals <- Key !!!
    inputsound = zeros(1,floor(dur/tdres)); 
    inputsound(1:2) = 1.;  % put 100 usec click at beginning of stimulus
    inputsound = aa * inputsound; 
    
    % Call up the auditory nerve model
    sout=anmod3m(CF, inputsound);
    
    subplot(5,1,iplot)
    plot((1:length(sout))*tdres*1e3,sout)
    set(gca,'Fontsize',14)
    xlabel('Time (msec)')
    ylabel('Sp/sec')
    text(15,max(sout)/2,['Peak SPL=' num2str(spl) 'dB SPL']);
    if iplot == 1
        title('Click Responses')
    end
end