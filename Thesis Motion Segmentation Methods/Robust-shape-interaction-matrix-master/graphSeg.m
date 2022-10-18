function graphSeg(calculatedSeg,groundTruthLabels)

    % Plot the ground truth labels
    subplot(2,1,1);
    gt_t = groundTruthLabels';
    plot(gt_t','.g');
    title('Ground Truth Labels');
    
    
    % Plot the calculated labels
    subplot(2,1,2);
    cs_t = calculatedSeg';
    plot(cs_t,'.b');
    
    % Plot the errors as red over the calculated labels
    hold on;
    errorIndexes = find(gt_t(:) ~= cs_t(:));
    plot(errorIndexes(:),cs_t(errorIndexes(:)),'.r');
    
    % Add some extra features for all graphs
    grid on;
    
end
