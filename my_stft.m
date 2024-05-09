function [spec, t, f] = my_stft(signal, fs, nch, hop_size)

    w_len = 256;
    nfft = 256;

    s_len = length(signal);
    n_frames = floor((s_len-w_len)/hop_size);
    
    y = zeros(nfft, n_frames, nch);
    for i = 1:n_frames
        i_samples = (i-1)*hop_size+1:(i-1)*hop_size+w_len;
        yi = signal(:, i_samples);
        Yi = fft(yi, nfft, 2);
        y(:, i, :) = Yi.';
    end
    
    t = (0:1:n_frames-1)*hop_size/fs;
    f = 0:fs/(nfft):(fs-fs/nfft)/2;
    spec = y(1:nfft/2,:, :);