function [GPCA_missrates,RANSAC_missrates,LSA_missrates,GPCA_times,RANSAC_times,LSA_times,GPCA_cgt,RANSAC_cgt,LSA_cgt] = multiview_multibody_affine_demo(path,sequencename,dataset_folder,GPCA_missrates,RANSAC_missrates,LSA_missrates,GPCA_times,RANSAC_times,LSA_times,GPCA_cgt,RANSAC_cgt,LSA_cgt)
% Lauch the included algorithms (GPCA with spectral clustering,
% LSA, RANSAC) on one sequence of the Hopkins 155 database
% 
% Inputs:
%   path            relative path from the current directory to
%                   the Hopkins 155 database
%   sequencename    name of the sequence
%
% Example
%   multiview_multibody_affine_demo('../Hopkins155','1R2RCR')

%move in the project directory
fullpath=fullfile(path,sequencename);
if(~exist(fullpath,'dir'))
    error(['Project directory ''' sequencename '''doesn''t exist'])
end

cdold=cd;
cd(fullpath)

% Alter data names
alter_data

X = reshape(permute(x(1:3,:,:),[1 3 2]),3*size(x,3),size(x,2));

%generate the data for the tests
generate_test_data

tic

% fprintf('GPCA with spectral clustering');
ts=cputime;
group = multiview_multibody_affine_spectral(yord,ngroups);
t=(cputime-ts)/size(group,1);

missrate=missclass(group,N,ngroups)/sum(N);
GPCA_missrates(end+1) = missrate;
% fprintf('\tMissclassification: %f%% Time elapsed: %fs\n\n',missrate*100,t);
GPCA_cgt{end+1} = group;

GPCA_times(end+1) = toc;

% graphSeg(group,s);
% title('GPCA Plots');
figure;
subplot(1,2,1); gscatter(X(1,:),X(2,:),s); axis equal; title('GroundTruth'); legend off; set(gca,'Color','k');
subplot(1,2,2); gscatter(X(1,:),X(2,:),group); axis equal; title('GPCA'); legend off; set(gca,'Color','k');


tic

% fprintf('RANSAC');
ts=cputime;
[b,group] = multiview_multibody_affine_ransac(yord,ngroups,0.00002);
t=(cputime-ts)/size(group,1);

missrate=missclass(group,N,ngroups)/sum(N);
RANSAC_missrates(end+1) = missrate;
% fprintf('\tMissclassification: %f%% Time elapsed: %fs\n\n',missrate*100,t);
RANSAC_cgt{end+1} = group;

RANSAC_times(end+1) = toc;

% figure;
% graphSeg(group,s);
% title('RANSAC Plots');

figure;
subplot(1,2,1); gscatter(X(1,:),X(2,:),s); axis equal; title('GroundTruth'); legend off; set(gca,'Color','k');
subplot(1,2,2); gscatter(X(1,:),X(2,:),group); axis equal; title('RANSAC'); legend off; set(gca,'Color','k');

tic

% fprintf('LSA');
ts=cputime;
group = multiview_multibody_affine_lsa(yord,ngroups,4*ngroups);
t=(cputime-ts)/size(group,1);

missrate=missclass(group,N,ngroups)/sum(N);
LSA_missrates(end+1) = missrate;
% fprintf('\tMissclassification: %f%% Time elapsed: %fs\n\n',missrate*100,t);
LSA_cgt{end+1} = group;

LSA_times(end+1) = toc;

% figure;
% graphSeg(group,s);
% title('LSA Plots');

figure;
subplot(1,2,1); gscatter(X(1,:),X(2,:),s); axis equal; title('GroundTruth'); legend off; set(gca,'Color','k');
subplot(1,2,2); gscatter(X(1,:),X(2,:),group); axis equal; title('LSA'); legend off; set(gca,'Color','k');

        
cd(cdold)
