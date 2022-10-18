% Moves a file into a newly created folder by the same name
clear; clc; 

addpath(pwd);

myDir = '/Users/christos/Desktop/Thesis Datasets/KT3DMoSeg'; %gets directory

cd(myDir); 

myFiles = dir(fullfile(myDir,'*.mat')); %gets all mat files in struct
for k = 1:length(myFiles)
    filename = myFiles(k).name;
    filenameWO_mat = filename(1:end-4);
    % fprintf(1, 'Now reading %s\n', filenameWO_mat);
    
    endStr = filenameWO_mat(end-5:end);
    % disp(endStr);
    
    if strcmp(endStr,'Tracks')
        newFolderName = filenameWO_mat;
        mkdir(newFolderName);
        movefile(filename,newFolderName);
    end
end