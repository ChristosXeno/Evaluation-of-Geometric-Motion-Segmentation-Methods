% Start from a blank workspace and screen
clear;
clc;
close all;
warning off;
% Define the file pathways taken by this method
method_folder_path = pwd;
datasets_folder = '/Users/christos/Desktop/Thesis Datasets';
% dataset_folder = '/Hopkins155'; % this can change
dataset_folder = '/Hopkins12'; 
% dataset_folder = '/MTPV62';
% dataset_folder = '/KT3DMoSeg';

% Aquire the dataset folder path
dataset_folder_path = strcat(datasets_folder,dataset_folder); % this can change

% Add the method's file to the path
addpath(method_folder_path);

% Load up a sequence from the dataset
file_path = dataset_folder_path;
cd(file_path);

% Select the current file with the dataset
file = dir;

% Create an array that stores the sequence names
sequence_names = {};

% Create an array that stores the calculated 
% ground truth values for each method used
GPCA_cgt = {};
RANSAC_cgt = {};
LSA_cgt = {};

% Create an array that stores the missrates for each method used
GPCA_missrates = [];
RANSAC_missrates = [];
LSA_missrates = [];

% Create an array that stores the times for each method used
GPCA_times = [];
RANSAC_times = [];
LSA_times = [];

% Loop through all sequences
for i = 1:length(file)
    if( (file(i).isdir == 1) && ~strcmp(file(i).name,'.') && ~strcmp(file(i).name,'..') )
		filepath = file(i).name;
		eval(['cd ' filepath]);
		
		f = dir;
        
        fprintf('\n%d/%d %s \n',i,length(file),filepath);
        
        for j = 1:length(f)
            %f(j).name
            if(~strcmp(f(j).name,'.') && ~strcmp(f(j).name,'..') && isempty(strfind(f(j).name,'.txt')) && isempty(strfind(f(j).name,'.png')) && isempty(strfind(f(j).name,'.jpg')))%(~isempty(strfind(f(j).name,'_truth.mat')) )			
                
                % Select the current sequence
                sequence_name = filepath;
                % sequence_name = '1R2TCR';
                % sequence_name = 'oc1R2RC_g12_Tracks';
                % sequence_name = 'PerspectiveNViewsData_flowershirt_Tracks';
                % sequence_name = 'Seq059_Clip01_Tracks';
                
                sequence_names{end+1} = sequence_name;

                % Run a demo for GPCA, RANSAC, LSA
                [GPCA_missrates,RANSAC_missrates,LSA_missrates,GPCA_times,RANSAC_times,LSA_times,GPCA_cgt,RANSAC_cgt,LSA_cgt] = multiview_multibody_affine_demo(file_path,sequence_name,dataset_folder,GPCA_missrates,RANSAC_missrates,LSA_missrates,GPCA_times,RANSAC_times,LSA_times,GPCA_cgt,RANSAC_cgt,LSA_cgt);
                
            end
        end
        % Move back out into the dataset folder
        cd ..
    end
end

% PerspectiveNViewsData_flowershirt_Tracks 

% Find the average and median missrates for the dataset
avgtol_GPCA = mean(GPCA_missrates);
avgtol_RANSAC = mean(RANSAC_missrates);
avgtol_LSA= mean(LSA_missrates);

medtol_GPCA = median(GPCA_missrates);
medtol_RANSAC = median(RANSAC_missrates);
medtol_LSA= median(LSA_missrates);


% Find min, max and avg times
mintime_GPCA = min(GPCA_times);
mintime_RANSAC = min(RANSAC_times);
mintime_LSA = min(LSA_times);

maxtime_GPCA = max(GPCA_times);
maxtime_RANSAC = max(RANSAC_times);
maxtime_LSA = max(LSA_times);

avgtime_GPCA = mean(GPCA_times);
avgtime_RANSAC = mean(RANSAC_times);
avgtime_LSA = mean(LSA_times);

% Display results
fprintf('\n');
disp(['Results on ' dataset_folder ':'])
disp(['GPCA:  Mean of all : ' num2str(100*avgtol_GPCA) '%' ', median of all: ' num2str(100*medtol_GPCA) '%;']);
disp(['RANSAC:  Mean of all: ' num2str(100*avgtol_RANSAC) '%' ', median of all: ' num2str(100*medtol_RANSAC) '%;']);
disp(['LSA:  Mean of all: ' num2str(100*avgtol_LSA) '%' ', median of all: ' num2str(100*medtol_LSA) '%;']);

disp(['GPCA:  Min Time : ' num2str(mintime_GPCA) ', Max Time: ' num2str(maxtime_GPCA) ', Avg Time: ' num2str(avgtime_GPCA)]);
disp(['RANSAC:  Min Time : ' num2str(mintime_RANSAC) ', Max Time: ' num2str(maxtime_RANSAC) ', Avg Time: ' num2str(avgtime_RANSAC)]);
disp(['LSA:  Min Time : ' num2str(mintime_LSA) ', Max Time: ' num2str(maxtime_LSA) ', Avg Time: ' num2str(avgtime_LSA)]);

% Move back to the folder where this method exists
cd(method_folder_path)

% Save Results
dataset_used = dataset_folder(2:end);

GPCA_results = {GPCA_missrates,GPCA_times,GPCA_cgt,sequence_names};
GPCA_variable_name = 'GPCA_results';
GPCA_save_mat_file = strcat(GPCA_variable_name,'_',dataset_used);

RANSAC_results = {RANSAC_missrates,RANSAC_times,RANSAC_cgt,sequence_names};
RANSAC_variable_name = 'RANSAC_results';
RANSAC_save_mat_file = strcat(RANSAC_variable_name,'_',dataset_used);

LSA_results = {LSA_missrates,LSA_times,LSA_cgt,sequence_names}; 
LSA_variable_name = 'LSA_results';
LSA_save_mat_file = strcat(LSA_variable_name,'_',dataset_used);

save(GPCA_save_mat_file,GPCA_variable_name);
save(RANSAC_save_mat_file,RANSAC_variable_name);
save(LSA_save_mat_file,LSA_variable_name);