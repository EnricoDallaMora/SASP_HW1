function [avg_pseudo_spec] = das_filter(y, fs, nch, theta, c, d)

    [spectrum, t_spec_axis, f_spec_axis] = my_stft(y, fs, nch);
    
    bands = f_spec_axis(4:4:end);
    
    a = zeros(nch, length(theta), length(bands));
    cov_est = zeros(nch, nch, length(bands));
    for f_c = bands
        a(:, :, bands==f_c) = exp(-1i*2*pi*f_c*d*sin(theta*pi/180).*(0:1:nch-1).'/c);
        for ii = 1:length(t_spec_axis)
            spec = squeeze(spectrum(f_spec_axis==f_c, ii, :));
            cov_est(:, :, bands==f_c) = cov_est(:, :, bands==f_c) + spec*spec'/length(t_spec_axis);
        end
    end
    
    p = zeros(length(bands), length(theta));
    for ii = 1:length(bands)
        for jj = 1:length(theta)
            p(ii, jj) = squeeze(a(:, jj, ii))'*cov_est(:, :, ii)*a(:, jj, ii)/nch^2;
        end
    end

    avg_pseudo_spec = sum(p, 1)/length(bands);