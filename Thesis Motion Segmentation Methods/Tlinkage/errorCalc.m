function error_rate = errorCalc(gt,calc)

    % Test if the segmentation starts with 0, if so
    % increase the value of all sementation methods
    if min(calc) == 0
        calc = calc + 1;
    end
    
    % Find the error rate
    error_rate = sum((gt(:)~=calc(:)))/length(gt);
    
end

