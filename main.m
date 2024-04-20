clc;
clear;
close all;

%% DATA AND INITIALIZATION
M = 16; % number of sensors
L = 0.45; % total length of array [m]
c = 343; % measured speed of sound [m/s]
d = L/(M - 1); % distance between sensors [m]

% d < lambda/2
lambda_min = 2*d;
f_max= c/lambda_min;  % anti-aliasing condition


[y, fs] = audioread('array_recordings.wav');

K = 512;
big_win = hann(K).';
big_hop = ceil(K*0.75);

N_frames = floor((length(y)-K)/big_hop);
theta = -90:1:90;
p_avg = zeros(N_frames, length(theta));
avg_theta = zeros(1, N_frames);

for kk = 1:N_frames
    y_w = y((kk-1)*big_hop+1:(kk-1)*big_hop + K, :).'.*big_win;
    
    % SIGNALS GENERATION
    
    w_len = 64; % length of window
    w = ones(1, w_len);
    R = 1;
    nfft = 512;
    
    
    % t = 0:1/fs:(length(y_st)-1)/fs;
    % limits = [0.9*min(mink(y_st, 1, 2)), 1.1*max(maxk(y_st, 1, 2))];
    % parfor i=1:M
    %     subplot(4, 4, i);
    %     plot(t, y_st(i, :));
    %     title(['Mic ', num2str(i)])
    %     ylim(limits)
    %     xlabel('time [s]')
    %     ylabel(['y_{', num2str(i), '}(t)'])
    %     grid on
    % end
    
    % PROCESSING
    [spectrum, t_spec_axis, f_spec_axis] = my_stft(y_w, fs, w, R, nfft, M);
    
    % [ff, tt] = meshgrid(t_spec_axis, f_spec_axis);
    % 
    % figure
    % parfor ii=1:M
    %     subplot(4, 4, ii)
    %     surf(ff, tt, abs(spectrum(:, :, ii)), 'EdgeAlpha', 0)
    %     view(2)
    %     title(['Spectrogram ', num2str(ii)])
    %     xlabel('time [s]')
    %     ylabel('frequency [Hz]')
    % end
    
    % Propagation vectors and covariance matrix estimate
    bands = f_spec_axis(4:4:end);
    a = zeros(M, length(theta), length(bands));
    cov_est = zeros(M, M, length(bands));
    
    for f_c = bands
        a(:, :, bands==f_c) = exp(-1i*2*pi*f_c*d*sin(theta*pi/180).*(0:1:M-1).'/c);
        for ii = 1:length(t_spec_axis)
            spec = squeeze(spectrum(f_spec_axis==f_c, ii, :));
            cov_est(:, :, bands==f_c) = cov_est(:, :, bands==f_c) + spec*spec'/length(t_spec_axis);
        end
    end
    
    % Pseudo-spectrum
    p = zeros(length(bands), length(theta));
        
    for ii = 1:length(bands)
        for jj = 1:length(theta)
            p(ii, jj) = squeeze(a(:, jj, ii))'*cov_est(:, :, ii)*a(:, jj, ii)/M^2;
        end
    end
    
    
    p_avg(kk, :) = sum(p, 1)/length(bands);
    
    [~, indexes] = max(abs(p_avg(kk, :)));
    
    avg_theta(kk) = theta(indexes);
    
    % close all
    % figure
    % for i = 1:length(bands)
    %     polarplot(theta*pi/180, abs(p(i, :)));
    %     title("Estimated DOA at frequency: " + bands(i) + " Hz");
    %     drawnow update;
    %     pause(0.1);
    % end
    
    % polarplot(theta*pi/180, abs(p_avg));
    % title("Average DOA across all frequency bands");
    % drawnow update;
    % pause(0.1);
end
%%
abs_p_avg = abs(p_avg);
abs_p_avg = abs_p_avg ./ max(abs_p_avg, [], 2);
[tax, theax] = meshgrid(theta, (0:1:N_frames-1)*big_hop/fs);
figure
surf(theax, tax, abs_p_avg, 'EdgeAlpha', 0)
view(2)
yticks([-90 -60 -30 0 30 60 90])
xticks(0:2:14)
xlabel('time [s]')
ylabel('DOA [deg]')
xlim([-0.5 15])
title('Normalized pseudo-spectrum in time')

%%

[u, v] = pol2cart(avg_theta*pi/180+pi/2, 0.25);
figure
polar(0, 0.25)
hold on
quiver(zeros(1, length(avg_theta)), zeros(1, length(avg_theta)), u, v, 'LineStyle','-', 'Color', 'green', 'LineWidth', 1.2)
hold on
stem(linspace(-0.45/2, 0.45/2, 16), zeros(1, 16),  'LineStyle', 'none', 'Marker','diamond', 'MarkerFaceColor','red', 'MarkerSize', 10, 'Color','black', 'LineWidth',2);
xlim([-0.25 0.25])
ylim([-0.1 0.4])
axis equal
grid on
%%
myVideo = VideoWriter('doas','MPEG-4');
myVideo.FrameRate = 30;
open(myVideo)
figure(100)
for i = 1:length(avg_theta)
    polar(0, 0.3);
    hold on
    stem(linspace(-0.45/2, 0.45/2, 16), zeros(1, 16),  'LineStyle', 'none', ...
        'Marker','diamond', 'MarkerFaceColor','red', 'MarkerSize', 10, ...
        'Color','black', 'LineWidth',2);
    hold on
    quiver(0, 0, u(i), v(i), 'LineStyle','-', 'Color', 'green', 'LineWidth', 1.2)
    hold off
    xlim([-0.25 0.25])
    ylim([-0.1 0.4])
    axis equal
    title(['Estimated DOA at time ', num2str(round(big_hop*i/fs, 2)), ' [s]'])
    xlabel('array elements coordinates [m]')
    xticks([-22.5e-2 -11.25e-2 0 11.25e-2 22.5e-2])
    yticks([])
    frame = getframe(100);
    writeVideo(myVideo, frame)
    drawnow limitrate
end
close(myVideo)