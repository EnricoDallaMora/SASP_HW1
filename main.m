%% DATA AND INITIALIZATION

clc;
clear all;

%[y,Fs] = audioread(rec.wav);

M=16; %number of sensors
L=0.45; %total length of array [m]
c=343; %measured speed of sound [m/s]
d=L/(M-1); %distance between sensors [m]

%% ANTI-ALIASING CONDITION 
%d<lambda/2
lambda_min=2*d;
f_max=c/lambda_min;  %anti-aliasing condition

%% SIGNALS GENERATION
fs=8000;
t_max=0.3;
t=0:t_max/fs:t_max;
y=zeros(M, length(t));
f0=10;
w_len=33; %length of window
w=hanning(w_len)'; %actual window
R=(w_len-1)/2+1; %hop size in samples (50%)
nfft=512;
n_frames=floor((length(y(1, :))-w_len)/R);
spectrum=zeros(M, n_frames, nfft);
theta=-90:1:90;
a=zeros(M, length(theta));


figure(1)
for i=1:M
    y(i,:)=square(2*pi*f0*t+i*pi/M);
    subplot(8, 2, i);
    plot(t, y(i, :));
    [spectrum(i, :, :), t_spec_axis, f_spec_axis]=my_stft(y(i, :), fs, w, R, nfft);
    for freq=1:length(f_spec_axis)
        a(i, :)=exp(-1i*2*pi*f_spec_axis(freq)*d*sin(theta*pi/180)*i/c);
        for frame=1:n_frames
            covariance_R(:, :)=1/n_frames*spectrum(:, :, freq)'*spectrum(i, frame, freq);
        end
    end
end


