% change EEG file sampling rate to 164

load(['perceived_env_3_4Hz.mat']); % or load short wav with 16000 fs x 3

% get syllable peaks
t = 1:size(envelope_filt,1);
TF = islocalmax(envelope_filt);
peak_times = t(TF);
number_peaks = length(peak_times);

% import data
test_data = EEG.data;
start = EEG.event(:,3).latency; 
ending = start+(164*60); % 1 minute
%ending = EEG.event(:,4).latency;
test_data = test_data(:,start:ending);

% parameters
n_w = number_peaks; % number of syllables

% make array of windows 
w_time = 1900/(1000/fs_new); % -800 to 1100 ms epochs % keep same window size as Etchell 
windows = zeros(n_w, round(w_time), 32);   
before = round(131.2);
after = ceil(180.4);

% loop over all syllables
for i = 1:n_w
    if peak_times(i)>before
        if peak_times(i) < size(test_data,2) - after
            onset = peak_times(i) - before;
            offset = peak_times(i) + after;
            windows(i,:,:)= test_data(:, onset:offset-1)';
        end
    end
end

% plot time domain
figure
plot(test_data(24,peak_times(10)-before:peak_times(10)+after))
title("Raw Signal for one epoch at channel Cz")

% get average of left channels
chan_avg_left = mean(windows(:,:,1:16), 3); %one window, all times minus ISI, all channels
chan_avg_trials_left = mean(chan_avg_left, 1);

% plot left hemisphere average
figure
plot(chan_avg_trials_left)
title('EEG Signal Averaged Across Windows and Averaged Across LH Channels')
xlabel('timestamps (1 timestamp = 6 ms (164 Hz Fs))')
ylabel('microVolt')

% get average of right channels
chan_avg_right = mean(windows(:,:,17:31), 3); % one window, all times, all channels
chan_avg_trials_right = mean(chan_avg_right, 1);

% plot right hemisphere average
figure
plot(chan_avg_trials_right)
title('EEG Signal Averaged Across Windows and Averaged Across RH Channels')
xlabel('timestamps (1 timestamp = 2 ms (500 Hz Fs))')
ylabel('microVolt')

%% Time-frequency analysis
num_frex = 20;    
eegpower_all_seqs = zeros(n_w,num_frex,round(w_time));                  % store results for each window, 20 frex

for k = 1:n_w
    % Define wavelet parameters
    left_chans = chan_avg_left(k,:);                    
    right_chans = chan_avg_right(k,:);
    
    min_freq =  1;
    max_freq = 40;     
    num_frex = 20;     
    
    time = -1:1/EEG.srate:1;    % size of wavelet
    frex = logspace(log10(min_freq),log10(max_freq),num_frex); % log if lower frequencies are of interest (up to 40/50 Hz) linear if higher
    s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
    
    % Definte convolution parameters
    n_wavelet            = length(time);
    n_data               = round(w_time);                   
    n_convolution        = n_wavelet+n_data-1;
    n_conv_pow2          = pow2(nextpow2(n_convolution));
    half_of_wavelet_size = (n_wavelet-1)/2;  % should be two periods of center frequency (20 Hz --> 100 ms)
    
    % get FFT of data
    eegfft = fft(right_chans,n_conv_pow2);                                           % change hemisphere!!
    
    % initialize
    eegpower = zeros(num_frex,round(w_time)); % if seq_time corresponds to eeg.pnts here, it has to be multiplied above
    baseidx = [1, 400];
    
    for fi=1:num_frex
        
        wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
        % convolution
        complex = wavelet.*eegfft;
        eegconv = ifft(complex);
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        
        % Average power over trials (baseline transform)
        temppower = abs(eegconv).^2;    %% mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2); %% trials
        eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2):2)));    
    end

    eegpower_all_seqs(k,:,:)= eegpower;  
end

% plot time-frequency result
A = eegpower_all_seqs(1:200,:,:);
A = nanmean(A);
A = reshape(A, [20, 312]);
A_mean = mean(A);
A_rescaled = A - A_mean;

% plot time frequency % minus 1500 to remove time between sequences
figure
contourf(1:round(w_time), frex,A_rescaled, max_freq, 'linecolor','none') %seq_time 
yline(12,'w-','lineWidth',2)
yline(15,'w-','LineWidth',2)
set(gca,'clim',[-2 2])
xlabel("time (ms)",'FontSize',16)
ylabel("Frequency (Hz)",'FontSize',16)
title('CWNS (9 y/o)', 'FontWeight','bold','FontSize',20)            % adjust title
subtitle('Synchronize - Central-Parietal Sites','FontSize',20)
xlim([57 262])
%xticklabels({'-400','-200','0','200','400','600','800'})
xline([81.9 131.2 180.5 229.8 279.1])
colormap turbo
colorbar

%% Get percentage signal change
% rescale to get %signal change of beta band
betaband = A_rescaled(12:15,:);  %12:15

avgFirst = mean(betaband,1);
avg2norm = sum(avgFirst)/size(betaband,2);
intermediate = mean(betaband,1) - avg2norm;
norm_signal = intermediate./max(abs(intermediate));

figure
plot(intermediate,'LineWidth',3, "Color", 'b')
xlabel("time (ms)",'FontSize',16)
xlim([57 262])
ylabel("%Signal Change",'FontSize',16)
title('CWNS (8 y/o)', 'FontWeight','bold','FontSize',20)        % adjust title
subtitle('Synchronize - Left Hemisphere: Beta Band (12-15 Hz)','FontSize',20)
%xticklabels({'-400','-200','0','200','400','600','800'})
xline([81.9 131.2 180.5 229.8 279.1])

% save result
save("Sync_eegpower.mat","A_rescaled","A","betaband")
