function mav = movingAverage(data, windowSize)
    % movingAverage calculates the moving average of the input data
    % using a specified window size.
    %
    % Args:
    %   data (vector): Input data vector.
    %   windowSize (integer): The size of the moving window.
    %
    % Returns:
    %   mav (vector): The moving average of the data.

    if length(data) < windowSize
        error('Window size must be less than or equal to the length of the data.');
    end

    % Initialize the output moving average vector
    mav = zeros(size(data));

    % Calculate the moving average using a sliding window
    for i = 1:length(data) - windowSize + 1
        mav(i) = mean(data(i:i + windowSize - 1));
    end

    % Handle the tail of the data where window does not fit
    for i = length(data) - windowSize + 2:length(data)
        mav(i) = mean(data(i:end));
    end

    % Alternatively, MATLAB's built-in 'smoothdata' function can be used
    % mav = smoothdata(data, 'movmean', windowSize);
end
