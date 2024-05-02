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

time = 0:1/fs:(length(y)-1)/fs;

figure
sgtitle('Array signals')
for ii = 1:16
    subplot(4, 4, ii)
    plot(time, y(:, ii)./max(maxk(y, 1, 2)))
    ylim([-1.1 1.1])
    xlabel('time [s]')
    xticks([2 6 10 14])
    title(['Mic ' num2str(ii)])
    grid on
end

%% PROCESSING
K = 1024;
big_win = hann(K).';
big_hop = floor(K*0.75);
N_frames = floor((length(y)-K)/big_hop);

theta = -90:1:90;
p_avg = zeros(N_frames, length(theta));
avg_theta = zeros(1, N_frames);

for kk = 1:N_frames
    y_w = y((kk-1)*big_hop+1:(kk-1)*big_hop + K, :).'.*big_win;    
    
    p_avg(kk, :) = das_filter(y_w, fs, M, theta, c, d);
    
    [~, indexes] = max(abs(p_avg(kk, :)));
    
    avg_theta(kk) = theta(indexes);
end

%% PSEUDO-SPECTRUM OVER TIME
abs_p_avg = abs(p_avg);
abs_p_avg = abs_p_avg ./ max(abs_p_avg, [], 2);
[tax, theax] = meshgrid(theta, (0:1:N_frames-1)*big_hop/fs);
figure
subplot(2, 1, 1)
surf(theax, tax, abs_p_avg, 'EdgeAlpha', 0)
view(2)
yticks([-90 -60 -30 0 30 60 90])
ylim([-90 90])
xticks(0:2:14)
xlim([0 14.6])
xlabel('time [s]')
ylabel('\theta [deg]')
title('Normalized pseudo-spectrum over time')
colorbar

subplot(2, 1, 2)
plot((0:1:N_frames-1)*big_hop/fs, avg_theta, 'LineWidth', 1.2)
xticks(0:2:14)
xlim([0 14.6])
yticks([-90 -60 -30 0 30 60 90])
ylim([-90 90])
xlabel('time [s]')
ylabel('\theta_{dir} [deg]')
title('Estimated DOA')
grid on

%% STATIC PLOT
[u, v] = pol2cart(avg_theta*pi/180+pi/2, 0.25);

x_c = -0.25:0.001:0.25;
circle = sqrt(-x_c.^2+(0.25)^2);

degrees = {'-90°', '-60°', '-30°', '0°', '30°', '60°', '90°'};

figure
plot(x_c, circle, 'Color', 'k', 'LineWidth', 1.1, 'LineStyle', '--')
for ll = 0:6
    line([0 0.27*cos(pi/6*ll)], [0 0.27*sin(pi/6*ll)], 'Color', 'k', 'LineWidth', 1.1, 'LineStyle', '--')
    text(0.29*cos(pi/6*ll), 0.29*sin(pi/6*ll), degrees(ll+1), 'HorizontalAlignment', 'center')
end
hold on
plot(linspace(-0.45/2, 0.45/2, 16), zeros(1, 16),  'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', ...
    'MarkerSize', 8);
hold on
quiver(zeros(1, length(avg_theta)), zeros(1, length(avg_theta)), u, v, 'LineStyle','-', ...
    'Color', 'r', 'LineWidth', 1.2)

xlim([-0.3 0.3])
ylim([-0.1 0.4])
yticks([])
xticks([-0.225 -0.1125 0 0.1125 0.225])
xlabel('array elements coordinates [m]')
title('Estimated DOAs for each time frame')
axis equal
grid minor

%% DYNAMIC PLOT
myVideo = VideoWriter('doas','MPEG-4');
myVideo.FrameRate = 20;
open(myVideo)
figure(100)
for i = 1:length(avg_theta)
    plot(x_c, circle, 'Color', 'k', 'LineWidth', 1.2, 'LineStyle', '--')
    for ll = 0:6
        line([0 0.27*cos(pi/6*ll)], [0 0.27*sin(pi/6*ll)], 'Color', 'k', 'LineWidth', 1.2, 'LineStyle', '--')
        text(0.29*cos(pi/6*ll), 0.29*sin(pi/6*ll), degrees(ll+1), 'HorizontalAlignment', 'center')
    end
    hold on
    plot(linspace(-0.45/2, 0.45/2, 16), zeros(1, 16),  'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 8)
    hold on
    quiver(0, 0, u(i), v(i), 'LineStyle','-', 'Color', 'red', 'LineWidth', 2, 'MaxHeadSize', 0.3)
    hold off
    xlim([-0.25 0.25])
    ylim([-0.1 0.4])
    axis equal
    title(['Estimated DOA at time ', num2str(big_hop*i/fs, '%.2f'), ' [s]'])
    xlabel('array elements coordinates [m]')
    xticks([-22.5e-2 -11.25e-2 0 11.25e-2 22.5e-2])
    yticks([])
    grid minor
    frame = getframe(100);
    writeVideo(myVideo, frame)
    drawnow limitrate
end
close(myVideo)