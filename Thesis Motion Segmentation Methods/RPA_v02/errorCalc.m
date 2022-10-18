function error_rate = errorCalc(gt,calc)

    % Find the error rate (SEEMS INCORECT CAUSE STUFF IS MISORDERED..)
    error_rate = sum((gt(:)~=calc(:)))/length(gt);
    
end

