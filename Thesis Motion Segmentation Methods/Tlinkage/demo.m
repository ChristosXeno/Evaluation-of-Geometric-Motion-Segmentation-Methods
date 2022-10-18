clear; clc; close all;
fprintf('Running T-Linkage: \n');

addpath(genpath('.'));
% This code just simply run the T-Linkage algorithm on the example data set
% "star5".
% Loading data: X contains data points, whereas G is the ground truth
% segmentation
% load 'star5.mat'; N = size(X,2);
% In order to work with a specific model, T-Linkage needs to be given:
% - distFun: distance between points and models
% - hpFun: returns an estimate model given cardmss points
% - fit_model: least square fitting function
% 
% In this example we want to estimate lines so distFun is the euclidean
% distance between a point from a line in the plane and cardmss=2.
% Other  possible models are 'line', 'circle',
% fundamental matrices ('fundamental') and 'subspace4' (look in 'model_spec' folder).
%



% Define the file pathways taken by this method
method_folder_path = pwd;
datasets_folder = '/Users/christos/Desktop/Thesis Datasets';
% dataset_folder = '/Hopkins155'; 
% dataset_folder = '/Hopkins12'; 
% dataset_folder = '/MTPV62';
dataset_folder = '/KT3DMoSeg';

% Aquire the dataset folder path
dataset_folder_path = strcat(datasets_folder,dataset_folder); % this can change

% Add the method's file to the path
addpath(method_folder_path);

% Load up a sequence from the dataset
file_path = dataset_folder_path;
cd(file_path);

% Select the current file with the dataset
file = dir;

% Create an array to store the error rates
error_rates = [];

% Create an array to store times
times = [];

% Create an array that stores the sequence names
sequence_names = {};

% Create an array that stores the calculated 
% ground truth values for each method used
cgt = {};

% Loop through all sequences
for i = 1:length(file)
    if( (file(i).isdir == 1) && ~strcmp(file(i).name,'.') && ~strcmp(file(i).name,'..') )
		filepath = file(i).name;
		eval(['cd ' filepath]);
		
		f = dir;
        
        fprintf('\n%d/%d %s \n',i,length(file),filepath);
        
        for j = 1:length(f)
            % f(j).name
            if(~strcmp(f(j).name,'.') && ~strcmp(f(j).name,'..') && isempty(strfind(f(j).name,'.txt')) && isempty(strfind(f(j).name,'.png')) && isempty(strfind(f(j).name,'.jpg')) )
                
                % Select the current sequence
                sequencename = filepath;
                % sequencename = 'cars3';
                % sequencename = 'oc1R2RC_g12_Tracks';
                % sequencename = 'MTPV_data_Boat_Tracks';
                % sequencename = 'Seq009_Clip03_Tracks';

                sequence_names{end+1} = sequencename;
                
                % Move in the project directory
                fullpath=fullfile(file_path,sequencename);
                if(~exist(fullpath,'dir'))
                    error(['Project directory ''' sequencename '''doesn''t exist'])
                end

                cdold=cd;
                cd(fullpath)

                % Alter data names
                alter_data
                
                % Separate and order data to find error rates
                generate_test_data
                
                X = reshape(permute(x(1:2,:,:),[1 3 2]),2*size(x,3),size(x,2));

                % X = reshape(permute(x(1:3,:,:),[1 3 2]),3*size(x,3),size(x,2));
                

                %%
                N = size(X,2);
                
                % Begin the timer
                tic

                [distFun, hpFun, fit_model, cardmss] = set_model('subspace4');
                %% Conceptual representation of points

                %T-linkage starts, as Ransac with random sampling:
                % Unform sampling can be adopted
                S = mssUniform(X, 5*N,cardmss); 
                % in order to reduce the number of hypotheses also a localized sampling can
                % be used:
                %
                %       D = pdist(X','euclidean');  D = squareform(D);
                %       S = mssNorm( X, D, 2*N, cardmss);
                %

                H = hpFun( X, S ); 
                % flat = fit_flat(X,3*size(X,3));

                % generating a pool of putative hypotheses H.
                % The residuals R between points and model
                R = res( X, H, distFun );
                % R = res_flat(X,flat);
                % are used for representing points in a conceptual space.
                % In particular a preference matrix P is built depicting by rows points
                % preferences.
                % 
                       % 1.3e-1, 4.5e-1, 4.3e-1
                epsilon= 1.3e-2; %An inlier threshold value  epsilon has to be specified.
                P  = prefMat( R, epsilon, 1 );

                %% Clustering

                %T-Linkage clustering follow a bottom up scheme in the preference space

                C = tlnk(P); % C = tlnk_fast(P);

                % C is a vector of labels, points belonging to the same models share the
                % same label.
                %% Outlier rejection step

                %T-Linkage fit a model to all the data points. Outlier can be found in
                %different ways (T-Linkage is agonostic about the outlier rejection strategy),
                %for example discarding too small cluster, or exploiting the randomness of
                %a model.
                C  = outlier_rejection_card( C, cardmss );
                
                cgt{end+1} = C;
                
                % C = res_flat(C', fit_flat(C', 0));
                % C = C';
                
                % Outliers are labelled by '0'

                %% Showing results
                figure
                subplot(1,2,1); gscatter(X(1,:),X(2,:),G); axis equal; title('Ground Truth'); legend off; set(gca,'Color','k');
                subplot(1,2,2); gscatter(X(1,:),X(2,:),C); axis equal; title('T-Linkage'); legend off; set(gca,'Color','k');

                % Find the error rate
                % error_rate = errorCalc(G,C);
                % error_rate = missclass(C,PPG,ngroups,G)/sum(PPG);
                error_rate = missclass(C,PPG,ngroups,G);
                
                % fprintf('Error Rate: %f\n',error_rate);
                error_rates(end+1) = error_rate;
                
                times(end+1) = toc;
            end
        end
        % Move back out into the dataset folder
        cd ..
    end
end

% Find the mean and median error rates
avgtol = mean(error_rates);
medtol = median(error_rates);

% Find the min, max and average (avg can be used to find total) times
mintime = min(times);
maxtime = max(times);
avgtime = mean(times);

fprintf('\n');
disp(['Results on ' dataset_folder ':'])
disp(['T-Linkage:  Mean of all : ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
disp(['Min Time : ' num2str(mintime) ', Max Time: ' num2str(maxtime) ', Avg Time: ' num2str(avgtime)]);

% Move back to the original folder
cd(method_folder_path);

% Save Results
dataset_used = dataset_folder(2:end);

TLinkage_results = {error_rates,times,cgt,sequence_names};
TLinkage_variable_name = 'TLinkage_results';
TLinkage_save_mat_file = strcat(TLinkage_variable_name,'_',dataset_used);

save(TLinkage_save_mat_file,TLinkage_variable_name);

%% Reference
% When using the code in your research work, please cite the following paper:
% Luca Magri, Andrea Fusiello, T-Linkage: A Continuous Relaxation of
% J-Linkage for Multi-Model Fitting, CVPR, 2014.
%
% For any comments, questions or suggestions about the code please contact
% luca (dot) magri (at) unimi (dot) it

