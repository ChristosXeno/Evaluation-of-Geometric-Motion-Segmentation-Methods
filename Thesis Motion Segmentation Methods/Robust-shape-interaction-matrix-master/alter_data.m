% Changes data to something similar to hopkins

% Test what dataset is being used
if strcmp(dataset_folder,'/Hopkins12')
    
    % Load the data structure
    eval(['load ' f(ind).name]);
    
    % Convert the data appropriately
    s = Data.GtLabel;                   % Ground truth
    frames = Data.nFrames;              % Frames
    points = Data.nSparsePoints;        % Points
    y = Data.ySparse;                   % y
    
    y_r1 = y(1,:,:);
    y_r2 = y(2,:,:);
    x_r1 = 2.*(y_r1(:,:,:)-min(y_r1))./(max(y_r1)-min(y_r1)) - 1;
    x_r2 = 2.*(y_r2(:,:,:)-min(y_r2))./(max(y_r2)-min(y_r2)) - 1;
    x = ones(3,length(y(1,:,1)),frames);      
    x(1,:,:) = x_r1;
    x(2,:,:) = x_r2;                    % x
    
elseif strcmp(dataset_folder,'/MTPV62')
    
    % Load the data structure
    eval(['load ' f(ind).name]);
    
    % Convert the data appropriately
    s = Data.GtLabel;               % Ground truth
    frames = Data.nFrames;              % Frames
    points = Data.nSparsePoints;        % Points
    y = Data.ySparse;                   % y
    
    y_r1 = y(1,:,:);
    y_r2 = y(2,:,:);
    x_r1 = 2.*(y_r1(:,:,:)-min(y_r1))./(max(y_r1)-min(y_r1)) - 1;
    x_r2 = 2.*(y_r2(:,:,:)-min(y_r2))./(max(y_r2)-min(y_r2)) - 1;
    x = ones(3,length(y(1,:,1)),frames);      
    x(1,:,:) = x_r1;
    x(2,:,:) = x_r2;                    % x
        
elseif strcmp(dataset_folder,'/KT3DMoSeg')
    
    % Load the data structure
    eval(['load ' f(ind).name]);
    
    % Convert the data appropriately
    s = Data.GtLabel;                   % Ground truth
    frames = Data.nFrames;              % Frames
    points = Data.nSparsePoints;        % Points
    y = Data.ySparse;                   % y
    
    y_r1 = y(1,:,:);
    y_r2 = y(2,:,:);
    x_r1 = 2.*(y_r1(:,:,:)-min(y_r1))./(max(y_r1)-min(y_r1)) - 1;
    x_r2 = 2.*(y_r2(:,:,:)-min(y_r2))./(max(y_r2)-min(y_r2)) - 1;
    x = ones(3,length(y(1,:,1)),frames);      
    x(1,:,:) = x_r1;
    x(2,:,:) = x_r2;                    % x
    
else 
    
    %load ground - truth
    eval(['load ' f(ind).name]);
    y_test = y(1,:,:);
end