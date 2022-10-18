clear; clc; close all;

addpath(genpath('.'));
warning off;
% Define the file pathways taken by this method
method_folder_path = pwd;
datasets_folder = '/Users/christos/Desktop/Thesis Datasets';
% dataset_folder = '/Hopkins155'; % this can change
% dataset_folder = '/Hopkins12'; 
dataset_folder = '/MTPV62';
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
            %f(j).name
            if(~strcmp(f(j).name,'.') && ~strcmp(f(j).name,'..') && isempty(strfind(f(j).name,'.txt')) && isempty(strfind(f(j).name,'.png')) && isempty(strfind(f(j).name,'.jpg')))%(~isempty(strfind(f(j).name,'_truth.mat')) )			
                
                % Select the current sequence
                sequencename = filepath;
                % sequencename = '1R2TCR';
                % sequencename = 'oc1R2RC_g12_Tracks';
                % sequencename = 'Hopkins50_two_cranes_Tracks';
                % sequencename = 'Seq009_Clip03_Tracks';

                sequence_names{end+1} = sequencename;

                % Move in the project directory
                fullpath=fullfile(file_path,sequencename);
                if(~exist(fullpath,'dir'))
                    error(['Project directory ''' sequencename '''doesn''t exist'])
                end

                cdold=cd;
                cd(fullpath)
                
                
                % Begin the timer
                tic
                
                % Begin the segmentation process
                % nameFM  = 'biscuitbookbox';
                % 
                model = 'fundamental';

                [ distFun, hpFun, fit_model, cardmss, idegen, d] = set_model( model );

                alter_data;
                generate_test_data;
                
                X = reshape(permute(y(1:2,:,:),[1 3 2]),2*size(y,3),size(y,2));

                % sequencename = nameFM;
                % out = 1;
                % [y,G,img] = load_data_FM(sequencename,out);
                % 
                % % normalization
                % [dat_img_1 T1] = normalise2dpts(y(1:3,:));
                % [dat_img_2 T2] = normalise2dpts(y(4:6,:));
                % X = [ dat_img_1 ; dat_img_2 ];


                N = numel(G);
                k = max(G); % number of models


                sigma_gt = 0.0054;


                %% Preference Trick


                % guided sampling

                w = 0.5;
                blk = 3*N;
                S  = mssWeighted( X, 6*N, blk, 'cauchy', model, w, sigma_gt);
                S=S(blk+1:end,:);



                H = hpFun(X,S); %hypotheses

                R = res( X, H, distFun ); disp('Residuals computed')

                P = prefMat(R, sigma_gt, 6); % preference matrix

                K = exp(- (squareform(pdist(P,@tanimoto))).^2);  % similarity matrix

                %% Robust PCA

                lambda = 1/sqrt(size(K,1));
                [K_rpca, E_hat, ~] = inexact_alm_rpca(K, lambda);   


                %% symmetric matrix factorization

                [Uinit, mekmeans]  = guess_init( K_rpca, k , G);


                params.Hinit=Uinit; params.maxiter = 100000;
                [U, iter, obj] = symnmf_anls(K_rpca, k,  params);
                indU = indMax(U);


                % segmentations obtained from snmf

                F = seg_from_binaryU(U);

                softIndU= U;    softIndU(indU==0)=0; % mlsac like


                %% Model extraction

                niter = 100;
                
                % niter = 1000;
                % niter = 10;

                [Phi, Z] =  rinforzino(X, S, P, F, softIndU , model, sigma_gt, niter);

                [~,I] = max(Phi'*softIndU,[],1);
                mss = Z(I,:);


                %% refinement using robust statistic

                cost = 1.5271;

                C = segmentation( mss, X, model, U , sigma_gt, cost ,'nearest');
                
                cgt{end+1} = C;

                %% visualization
                figure;
                subplot(1,2,1); gscatter(X(1,:),X(2,:),G); axis equal; title('GroundTruth'); 
                subplot(1,2,2); gscatter(X(1,:),X(2,:),C); axis equal; title('RPA'); 

                % Find the error rate
                % error_rate = errorCalc(G,C);
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
disp(['RPA:  Mean of all : ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
disp(['Min Time : ' num2str(mintime) ', Max Time: ' num2str(maxtime) ', Avg Time: ' num2str(avgtime)]);

% Move back to the original folder
cd(method_folder_path);

% Save Results
dataset_used = dataset_folder(2:end);

RPA_results = {error_rates,times,cgt,sequence_names};
RPA_variable_name = 'RPA_results';
RPA_save_mat_file = strcat(RPA_variable_name,'_',dataset_used);

save(RPA_save_mat_file,RPA_variable_name);

