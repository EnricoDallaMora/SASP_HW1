function [spec, t, f]=my_stft(signal, fs, window, hop_size, nfft)

s_len=length(signal);
w_len=length(window);
n_frames=floor((s_len-w_len)/hop_size);

y=zeros(n_frames, nfft);
for i=1:n_frames
    i_samples=(i-1)*hop_size+1:(i-1)*hop_size+w_len;
    yi=signal(i_samples).*window;
    Yi=fft(yi, nfft);
    y(i, :)=y(i, :)+Yi; 
end

time=s_len/fs;
t=0:time/(n_frames-1):time;
f=0:fs/(nfft-1):fs;
spec=y;