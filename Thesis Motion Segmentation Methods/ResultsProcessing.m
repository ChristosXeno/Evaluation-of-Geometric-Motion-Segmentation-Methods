%%
% This processes the Results Aquired from all methods

% Start from a blank workspace and screen
clear; clc; close all;

% Define the file pathways taken by this method (change these if needed)
current_folder = '/Users/christos/Desktop/Thesis Motion Segmentation Methods';
datasets_folder = '/Users/christos/Desktop/Thesis Datasets';

% Add the current file to the path
addpath(current_folder);

% Determine the method
% method = 'GPCA';
% method = 'LSA';
method = 'RSIM';
% method = 'LRSSC';
% method = 'RANSAC';
% method = 'TLinkage';
% method = 'RPA';

% Define the dataset
dataset = 'Hopkins155';
% dataset = 'Hopkins12'; 
% dataset = 'MTPV62';
% dataset = 'KT3DMoSeg';

% Define the dataset folder
dataset_folder = strcat('/',dataset);

% Define the dataset folder path
dataset_folder_path = strcat(datasets_folder,'/',dataset);

% Aquire the results file name
results_struct_name = strcat(method,'_results');
results_file = strcat(results_struct_name,'_',dataset,'.mat');

% Load the results
results_struct = load(results_file);
method_results = results_struct.(results_struct_name);

% Name each cell of the results cell array
error_rates = method_results{1};
times = method_results{2};
cgts = method_results{3};
seq_names = method_results{4};


%% Across Dataset
% Find the mean and median error rates
avgerr = mean(error_rates);
mederr = median(error_rates);

% Find the min, max and average (avg can be used to find total) times
mintime = min(times);
maxtime = max(times);
avgtime = mean(times);

fprintf('\n');
disp(['Results on ' dataset_folder ':'])
disp([method ':  Mean of all : ' num2str(100*avgerr) '%' ', median of all: ' num2str(100*mederr) '%;']);
disp(['Min Time : ' num2str(mintime) ', Max Time: ' num2str(maxtime) ', Avg Time: ' num2str(avgtime)]);

%% Per Seq
% % Find a sequence according to it's name
seq_name_goal = '1R2RC';
seq_index = find(ismember(seq_names,seq_name_goal));
%seq_index = 1;

% seq_names(find(error_rates == max(error_rates)))
%seq_names(find(error_rates == min(error_rates)))

% Define the index of the sequence which is to be analysed (if needed)
seq_index_start = seq_index;

% Define the amount of sequences used
seqs = 1;
% seqs = length(error_rates);

% Define loop bounds 
if seqs == 1
    i_s = seq_index_start;
    i_e = seq_index_start; 
elseif (seqs == length(error_rates))
    i_s = 1;
    i_e = length(error_rates);
else
    i_s = seq_index_start;
    i_e = seq_index_start + seqs;
end

for seq = i_s:i_e

    % Aquire the error metric for one sequence
    cgt_seqs = cgts(seq);
    name_seqs = seq_names(seq);

    err_seq = error_rates(seq);
    time_seq = times(seq);
    cgt_seq = cgt_seqs{1};
    name_seq = name_seqs{1};


    file_path = dataset_folder_path;
    sequencename = name_seq;

    %%
    % Load up a sequence from the dataset
    cd(file_path);
    eval(['cd ' sequencename]);

    % Select the current file with the dataset
    file = dir;


    % Process the data
    fullpath=fullfile(file_path,sequencename);

    if(~exist(fullpath,'dir'))
        error(['Project directory ''' sequencename '''doesn''t exist'])
    end

    cdold=cd;
    cd(fullpath)

    alter_data

    X = reshape(permute(x(1:3,:,:),[1 3 2]),3*size(x,3),size(x,2));

    Y = reshape(permute(y(1:3,:,:),[1 3 2]),3*size(y,3),size(y,2));

    cd(current_folder);
    %%
    figure;
    fig_title = ['Results on: ', dataset_folder,': ', name_seq];
    sgtitle(fig_title,'Interpreter','none');
    subplot(1,2,1); gscatter(Y(1,:),Y(2,:),s); axis equal; title('GroundTruth'); legend off; set(gca,'Color','k');
    subplot(1,2,2); gscatter(Y(1,:),Y(2,:),cgt_seq); axis equal; title(method); legend off; set(gca,'Color','k');
    
    fprintf('\n');
    disp(['Error Rate on: ', dataset_folder,': ', name_seq, ': ',num2str(error_rates(seq_index)*100),'%']);
end

