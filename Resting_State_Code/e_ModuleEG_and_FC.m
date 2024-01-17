%% ESTMATE_NETWROK_METRICS 
% Estimation of brain network efficiency and modularity.                                                             
%                                                                         
% Requires Brain Connectivity Toolbox and functional connectivity maps.   
%                                                                         
% Julia Linke May/2020, edited April/2022                                               
%
%% Configuration
clear
clc

fprintf('   *** Estimation of module efficiency and connectivity ...\n')

bcn_dir = ''; %path to the brain connectivity toolbox (download here: https://sites.google.com/site/bctnet/)
data_dir = '.../netmats'; %path to the output of the previous script b_Netmats_216nodes
list_dir = '.../lists/';  %path to the list with the participant IDs
cocl_dir = '.../CoClass/'; %path to output from previous code (d_Coclass.m)
metric_dir = '.../Metrics';

addpath(bcn_dir)
addpath(data_dir)
addpath(list_dir)
addpath(metric_dir)

ID = fileread(strcat(list_dir,'GraphTheorySample.txt'));
ID = strsplit(ID);
ID(end)=[]

data = {'ses-1_task-rest'};

Res1EG = cell(length(data),1); for c = 1:numel(Res1EG), Res1EG{c} = cell(1,1);end
Res1Eloc = cell(length(data),1); for c = 1:numel(Res1Eloc), Res1Eloc{c} = cell(1,1);end
Res1With = cell(length(data),1); for c = 1:numel(Res1With), Res1With{c} = cell(1,1);end
Res1Btw = cell(length(data),1); for c = 1:numel(Res1Btw), Res1Btw{c} = cell(1,1);end

ii=0; 
jj=0;

parcellation_code   = 'Schaefer200';
no_of_nodes         = 216;
wlet_scale          = 5; 
FisherZ             = 1;
method_of_thr       = 'rnk';
res_par             = 1;


%% load Coclassification for 10% threshold
cd (cocl_dir)
restref=load('Withprior_ses-1_task-rest_scale5_rnk10_tau75_grouptau50.mat');

%indices for networks of interest
%might vary by rank-threshold, need to check and manually update
rest1 = {1, 2, 3, 4, 5, 8};

%% Main loop

   for s = 1:length(ID)

        for d = 1:length(data)

            disp(['       Subject: ' ID{s} ', Data: ' data{d}])

            % Load the connectivity matrix
            cd (fullfile(data_dir,ID{s}))
            W = load(strcat(ID{s}, '_' ,data{d} ,'_netmat-full_atlas-Schaefer2018-200P+17N_space-T1w.csv'));

            % Set the diagonal to zero
            W   = W - diag(diag(W));

            % Fisher's Z transform
            if FisherZ
                W = atanh(W); % Fisher's Z ranges between -inf and +inf
            end

            r_vec = W(find(triu(W,1)));
            RawFC = sum(r_vec)/length(r_vec);

            % Thresholding
            if strcmp(method_of_thr, 'rnk') % Rank thresholding
                r_up        = W(find(triu(W) & ~eye(no_of_nodes)));
                [Y, I]       = sort(abs(r_up), 'descend');
                index_thres = round(max(I)*0.1); %set it to 10 percent here
                r_thres     = Y(index_thres);

                W(abs(W)<=abs(r_thres)) = 0;

            elseif strcmp(method_of_thr, 'abs') % Absolute thresholding
                W(abs(W)<thr(k)) = 0; 
            end
            A = (ones(no_of_nodes, no_of_nodes).*(W ~= 0)); % Binarize
            
            %extract Eg and meanEloc for resting-state networks
             
            jj=0; 
            for  m = 1:length(rest1) 
                [index,]= find(restref.Ccm == rest1{m});
                A1=A(index,index);
                Eg      = efficiency_bin(A1);
                Eloc    = efficiency_bin(A1,1);
                MeanEloc = sum(Eloc)/no_of_nodes;
                Rest1EG{d,1}{s,m} = Eg;
                Rest1Eloc{d,1}{s,m} = MeanEloc;
                
                ii=0;
                for n = 1:length(rest1)
                    ii = ii+1;
                    if n == m
                           MM = W(find(restref.Ccm == rest1{m}), find(restref.Ccm == rest1{n}));
                           MM_up = MM(find(triu(MM) & ~eye(length(MM))));
                           Rest1With{d,1}{s,ii} = sum(MM_up)/(size(MM, 1)*size(MM, 2)); % According to Garcia et al. (2018)
                    elseif n>m
                           jj = jj+1;
                           MN = W(find(restref.Ccm == rest1{m}), find(restref.Ccm == rest1{n}));
                           MN_vec = reshape(MN, size(MN, 1)*size(MN, 2), 1); 
                           Rest1Btw{d,1}{s,jj} = sum(MN_vec)/(length(find(restref.Ccm == rest1{m}))*length(find(restref.Ccm == rest1{n})));
                    end
                end
                
            end
            
                       
       fprintf('       DONE.\n')

        end
    end


cd(cocl_dir)
%save('ModuleDetails10Perc.mat', 'Rest1EG', 'Rest1Eloc','Rest1With', 'Rest1Btw')
writecell(Rest1EG{1, 1}, 'ModuleDetails10Perc_Rest1EG.csv')
writecell(Rest1Eloc{1, 1}, 'ModuleDetails10Perc_Rest1Eloc.csv')
writecell(Rest1With{1, 1}, 'ModuleDetails10Perc_Rest1With.csv')
writecell(Rest1Btw{1, 1}, 'ModuleDetails10Perc_Rest1Btw.csv')


fprintf('   *** End of detail analysis. \n')
fprintf(['   ' datestr(now, 0) ' \n'])

%% End of the script