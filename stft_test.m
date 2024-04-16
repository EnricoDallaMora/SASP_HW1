fs=8000;
t_max=0.5;
t=0:1/fs:t_max;
f0=200;
f1=5;
x=square(2*pi*f0*t)+sin(2*pi*f1*t);
figure(1)
plot(t, x);
M=33; %length of window
w=hanning(M)'; %actual window
R=(M-1)/2+1; %hop size in samples (50%)
nfft=512;
[X, t, f]=my_stft(x, fs, w, R, nfft);
figure(2)
waterfall(f,t,real(X));
set(gca,'YDir', 'reverse')