% prepare data
% upsample EEG signal in EEGLab to match original audio
% Tools > Change sampling rate > 16000

% load data
load(['perceived_env_3_4Hz_short.mat']);
audio_all = envelope_filt'; 
audio_all = cat(1,audio_all, audio_all, audio_all, audio_all, audio_all, audio_all, audio_all); %cat(1, or 2?)

data = EEG.data; 
Fs = 16000;
onset = EEG.event(:,2).latency; 
w = 1; 
t = 140/w; 

all_fronts = zeros(1,w);
all_temps = zeros(1,w);

for i = 1:t
i
m = i-1;
onset = EEG.event(:,2).latency + m*(16000*w);
offset = onset + (16000*w)-1; 

% permutation 
r = randi([1 t*Fs],1,1000);

onset_perm = onset + r(v);
offset_perm = onset+ ((offset-onset)-(offset-onset_perm));

% set up parameters
window = round(Fs/10); %size of window (i.e. 779)
noverlap = round(Fs/20); %size of overlap (50% of window = 390)
nfft = 2048; %Fourier transform

% look at data    
audio = audio_all((Fs*w)*m+1:i*(w*Fs),1);
coh_vals = zeros(1,32);

for k = 1:32
    neural = data(k,onset:offset);
    freqFilt = [2.4, 4.4];
    neural = bandpass(neural,freqFilt,16000); % filter neural data
    
    % compute cross and auto-spectral density 
    [Pxx,f_x] = cpsd(audio, audio, window, noverlap, nfft, 10, 'centered'); %autospectra RMas %Fs?
    [Pyy,f_y] = cpsd(neural, neural, window, noverlap, nfft, 10, 'centered'); %autospectra LMas
    [Pxy,f_xy] = cpsd(audio, neural, window, noverlap, nfft, 10, 'centered'); %cross-spectra RMas x LMas
    
    % compute coherence 
    coh_xy = (abs(Pxy).^2)./(Pxx.*Pyy);
    coh = mean(coh_xy);
    coh_vals(1,k) = coh;
end

% frontal ROI
front_locs = [1,2,3,4,6,7,28,29,30,31,32]; 
front_vals = zeros(1,length(front_locs));

for k = 1:length(front_locs)
    front_vals(1,k) = coh_vals(1,front_locs(k)); 
end

frontal_coherence = mean(front_vals); 
all_fronts(1,i) = frontal_coherence;
%disp(frontal_coherence)

% temporal ROI
temp_locs = [15, 8, 14, 25, 19, 20]; 
temp_vals = zeros(1,length(temp_locs));

for k = 1:length(temp_locs)
    temp_vals(1,k) = coh_vals(1,temp_locs(k)); 
end

temporal_coherence = mean(temp_vals); 
all_temps(1,i) = temporal_coherence;
%disp(temporal_coherence)
end

front_coherence_avg = mean(all_fronts); 
temp_coherence_avg = mean(all_temps);

figure;
plot(all_fronts);
hold on
plot(all_temps);
legend('frontal', 'temporal')

disp(front_coherence_avg)
disp(temp_coherence_avg)

disp(size(neural))

% frontal ROI figure
plot_data = zeros(length(front_locs),length(b));

for k = 1:length(front_locs)
    plot_data(k,:) = data(k, onset:offset); 
end
plot_data = mean(plot_data);
plot_data =bandpass(plot_data,freqFilt,Fs); % filter neural data

figure;
plot(a);
hold on
plot(plot_data);
set(gca,'xlim',[30*164 33*164])
legend('Perceived Envelope','Frontal Entrainment')
