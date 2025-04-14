function Yt = get_Yt_window(OriSignal, winsize, morder)
% OriSignal: dimension × total_time
% winsize: number of time points per window
% morder: number of lags (VAR model order)
%
% Output:
% Yt: (dimension*(morder+1)) × num_windows × winsize

[dim, T] = size(OriSignal);

% First window starts at morder+1
start_index = morder + 1;

% Number of full windows that can be extracted
num_windows = floor((T - morder) / winsize);

% Preallocate output array
Yt = zeros(dim * (morder + 1), num_windows, winsize);

for w = 1:num_windows
    % Compute the time index of the first point in the current window
    win_start = start_index + (w - 1) * winsize;
    
    for t = 1:winsize
        current_t = win_start + t - 1;
        if current_t > T
            break; % Prevent overrun (shouldn't happen if num_windows is correct)
        end
        
        % Stack x(t), x(t-1), ..., x(t-morder)
        stacked = [];
        for lag = 0:morder
            stacked = [stacked; OriSignal(:, current_t - lag)];
        end
        
        Yt(:, w, t) = stacked;
    end
end

end
