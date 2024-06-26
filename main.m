%% DATA AND INITIALIZATION

clc;
clear;
close all;


[y, fs] = audioread('array_recordings.wav');
y = y.';
M = 16; %number of sensors
L = 0.45; %total length of array [m]
c = 343; %measured speed of sound [m/s]
d = L/(M - 1); %distance between sensors [m]

%% ANTI-ALIASING CONDITION 
%d<lambda/2
lambda_min = 2*d;
f_max= c/lambda_min;  %anti-aliasing condition

%% SIGNALS GENERATION
% 

% w = hanning(w_len)'; %actual window
% R = (w_len-1)/2 + 1; %hop size in samples (50%)
w_len = 128; %length of window
w = ones(1, w_len);
R = 1*w_len; %in samples, not in percentage
nfft = 512;
theta = -90:1:90;

% t_max = 0.3;
% t = 0:1/fs:t_max;
% f0 = 500;
% y = zeros(M, length(t));
% figure
% for i=1:M
%     y(i,:) = square(2*pi*f0*t+i*pi/(M));
%     subplot(4, 4, i);
%     plot(t, y(i, :));
%     title(['Mic ', num2str(i)])
%     ylim([-1.2 1.2])
% end

%% PROCESSING
[spectrum, t_spec_axis, f_spec_axis] = my_stft(y, fs, w, R, nfft, M);

[ff, tt] = meshgrid(t_spec_axis, f_spec_axis);

figure
for ii=1:M
    subplot(4, 4, ii)
    surf(ff, tt, abs(spectrum(:, :, ii)), 'EdgeAlpha', 0)
    view(2)
    title(['Spectrogram ', num2str(ii)])
end
%% Propagation vectors and covariance matrix estimate
bands = f_spec_axis(4:4:end);
temp_res=50;
a = zeros(M, length(theta), length(bands));
cov_est = zeros(M, M, length(bands), temp_res);
total_time_samples=length(t_spec_axis);
time_step=floor(total_time_samples/temp_res);

for f_c = bands
    a(:, :, bands==f_c) = exp(-1i*2*pi*f_c*d*sin(theta*pi/180).*(0:1:M-1).'/c);
    for tt = 1:temp_res
        for ii = 1+(tt-1)*time_step:tt*time_step
            spec = squeeze(spectrum(f_spec_axis==f_c, ii, :));
            cov_est(:, :, bands==f_c, tt) = cov_est(:, :, bands==f_c, tt) + spec*spec'/(total_time_samples/temp_res);
        end
    end
end

%% Pseudo-spectrum

p = zeros(length(bands), length(theta), temp_res);

for tt = 1:temp_res
    for ii = 1:length(bands)
        for jj = 1:length(theta)
            p(ii, jj, tt) = squeeze(a(:, jj, ii))'*cov_est(:, :, ii, tt)*a(:, jj, ii)/M^2;
        end
    end
end

p_avg = sum(p, 1)/length(bands);

[~, indexes] = max(abs(p_avg));

avg_theta = theta(indexes);


figure
for i = 1:temp_res
    polarplot(theta*pi/180, abs(p_avg(1, :, i)));
    title("Average DOA across all frequency bands, time: "+t_spec_axis(i*time_step)+" seconds");
    drawnow update;
    pause(0.1);
end

% figure
% for i = 1:length(bands)
%     polarplot(theta*pi/180, abs(p(i, :, 1)));
%     title("Estimated DOA at freqency: " + bands(i) + " Hz");
%     drawnow update;
%     pause(0.1);
% end