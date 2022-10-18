clear; clc; close all;

% Before running this code, SSC and LRR codes from the authors' websites
% should be in the subfolders with the following titles
% addpath SSC_ADMM_v1.1
% addpath code2
% addpath '/Users/christos/Desktop/greedysc_codes/supp_material_RSCT/include'
% addpath '/Users/christos/Desktop/greedysc_codes/SSC'
% addpath '/Users/christos/Desktop/greedysc_codes/LRR'

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

i_data = 0;

file = dir;

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
				
                alter_data
                
				if(max(s)==5)
					foundValidData = false;
				end
                break
            end
        end        
        cd ..
		
		if(foundValidData)

            % disp([filepath ': ' num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
            
            sequence_names{end+1} = filepath;
            
            fprintf('\n%d/%d %s \n',i,length(file),filepath);
            
            
            L = max(s);

            i_data = i_data + 1;

            N = size(x,2);
            F = size(x,3);
            p = 2*F;
            Y = reshape(permute(x(1:2,:,:),[1 3 2]),p,N);

            A0 = s; 
            [~,I] = sort(A0,1);

%                 %% SSC
%                 fprintf('\nRunning SSC..\n'); i_algo = 3;
%                 r = 0; affine = true; outlier = false; rho = 0.7; alpha = 800;
%                 tic
%                 [Z,A] = SSC(Y,r,affine,alpha,outlier,rho,s);
% 
%                 ET(i_algo,i_data)  = toc;
%                 CE(i_algo,i_data)  = computeCE(A,A0);

            %% LRSSC
            fprintf('Running LRSSC..\n');
            % LRR
            % The following scripts are copied from the LRR source code.
            tic
            lambda = 4;
            %run lrr
            Z = solve_lrr(Y,lambda);
            %post processing
            [U,S,V] = svd(Z,'econ');
            S = diag(S);
            r = sum(S>1e-4*S(1));
            U = U(:,1:r);S = S(1:r);
            U = U*diag(sqrt(S));
            %U = normr(U);
            U = U./repmat(sqrt(sum(U.^2,2)),1,size(U,2));
            LL = (U*U').^4;
            % spectral clustering
            D = diag(1./sqrt(sum(LL,2)));
            LL = D*LL*D;
            [U,S,V] = svd(LL);
            V = U(:,1:L);
            V = D*V;
            A = kmeans(V,L,'emptyaction','singleton','replicates',10,'display','off');

            % ET(i_algo,i_data)  = toc;
            % CE(i_algo,i_data)  = computeCE(A,A0);

            % SSC
            r = 0; affine = true; outlier = false; rho = 0.7; alpha = 800;
            % tic
            [Z,A] = SSC(Z,r,affine,alpha,outlier,rho,s);

            ET(i_data)  = toc;
            CE(i_data)  = computeCE(A,A0);
            
            cgt{end+1} = A;
            
            % figure;
            % graphSeg(A,A0);
            % title('LRSSC Labels');
            
            figure;
            subplot(1,2,1); gscatter(Y(1,:),Y(2,:),A0); axis equal; title('GroundTruth'); legend off; set(gca,'Color','k');
            subplot(1,2,2); gscatter(Y(1,:),Y(2,:),A); axis equal; title('LRSSC'); legend off; set(gca,'Color','k');

            
		end
	end
end

cd(method_folder_path);

meanCE = mean(CE')
medianCE = median(CE')
maxCE = max(CE')

meanET = mean(ET')
minET = min(ET')
maxET = max(ET')

%% Save Results
dataset_used = dataset_folder(2:end);

LRSSC_results = {CE,ET,cgt,sequence_names};
LRSSC_variable_name = 'LRSSC_results';
LRSSC_save_mat_file = strcat(LRSSC_variable_name,'_',dataset_used);

save(LRSSC_save_mat_file,LRSSC_variable_name);


