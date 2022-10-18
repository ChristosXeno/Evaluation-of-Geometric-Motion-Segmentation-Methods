% Motion segmentation with Robust Shape Interaction Matrix Method
% Pan Ji, pan.ji@anu.edu.au
% Nov 2014, @ANU
clear; close all; clc;
addpath(genpath(pwd))
method_folder_path = pwd;
datasets_folder = '/Users/christos/Desktop/Thesis Datasets';

dataset_folder = '/Hopkins155';
% dataset_folder = '/Hopkins12';
% dataset_folder = '/MTPV62';
% dataset_folder = '/KT3DMoSeg';

dataset_folder_path = strcat(datasets_folder,dataset_folder); % this can change
cd(dataset_folder_path) % cd to the dataset file path

warning off

file = dir;
ii = 0;
ii2 = 0;
ii3 = 0;
times = [];

% Create an array that stores the sequence names
sequence_names = {};

% Create an array that stores the calculated 
% ground truth values for each method used
cgt = {};

for i = 1:length(file)
	if( (file(i).isdir == 1) && ~strcmp(file(i).name,'.') && ~strcmp(file(i).name,'..') )
		filepath = file(i).name;
		eval(['cd ' filepath]);
		
		f = dir;
        foundValidData = false;
        for j = 1:length(f)
            %f(j).name
            if(~strcmp(f(j).name,'.') && ~strcmp(f(j).name,'..') && isempty(strfind(f(j).name,'.txt')) )%(~isempty(strfind(f(j).name,'_truth.mat')) )			
                ind = j;
                foundValidData = true;
                
                sequence_names{end+1} = filepath;
				
                alter_data
                
				if(max(s)==5)
					foundValidData = false;
				end
                break
            end
        end        
        cd ..
		
		if(foundValidData)
			N = size(x,2);
			F = size(x,3);
			D = 3*F;						
									
			X = reshape(permute(x(1:3,:,:),[1 3 2]),D,N);	% note here the all-one rows are also included
			
            tic
            
			[missrate, grp, bestRank, minNcutValue,W] = RSIM(X, s);	
            
            cgt{end+1} = grp;
            
            times(end+1) = toc;

			ii = ii+1;	
			Missrate(ii) = missrate;					
			disp([filepath ': ' num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
			if(max(s)==2)
				ii2 = ii2+1;
				Missrate2(ii2) = Missrate(ii);				
			else
				ii3 = ii3+1;
				Missrate3(ii3) = Missrate(ii);				
			end
		end
	end
end
% time = toc;
% avgtime = time/ii

% Find Missrates
avgtol = mean(Missrate);
medtol = median(Missrate);
avgtwo = mean(Missrate2);
medtwo = median(Missrate2);
avgthree = mean(Missrate3);
medthree = median(Missrate3);

% Find the min, max and average (avg can be used to find total) times
mintime = min(times);
maxtime = max(times);
avgtime = mean(times);

fprintf('\n');
disp(['Results on ' dataset_folder ':'])
disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%.']);
disp(['Min Time : ' num2str(mintime) ', Max Time: ' num2str(maxtime) ', Avg Time: ' num2str(avgtime)]);

% Move back to the original folder
cd(method_folder_path);

% Save Results
dataset_used = dataset_folder(2:end);

RSIM_results = {Missrate,times,cgt,sequence_names};
RSIM_variable_name = 'RSIM_results';
RSIM_save_mat_file = strcat(RSIM_variable_name,'_',dataset_used);

save(RSIM_save_mat_file,RSIM_variable_name);
