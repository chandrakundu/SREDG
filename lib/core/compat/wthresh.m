function y = wthresh(x, sorh, t)
%WTHRESH Soft or hard thresholding, the same as the Wavelet Toolbox function.

switch sorh
    case 'h'
        y = x .* (abs(x) > t);
    case 's'
        tmp = abs(x) - t;
        y = sign(x) .* ((tmp + abs(tmp)) / 2);
    otherwise
        error('wthresh:Type', 'Invalid threshold type ''%s''; use ''s'' or ''h''.', sorh);
end
end
