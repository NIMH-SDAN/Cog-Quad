%% ESTMATE_NETWROK_METRICS 
% Estimation of brain network efficiency and modularity.                                                             
%                                                                         
% Requires Brain Connectivity Toolbox and functional connectivity maps.   
%                                                                         
% Julia Linke April/2022                                               
%
%% Configuration
clear
clc

fprintf('   *** Estimation of network efficiency and modularity ...\n')

bcn_dir = ''; %path to the brain connectivity toolbox (download here: https://sites.google.com/site/bctnet/)
data_dir = '.../netmats'; %path to the output of the previous script b_Netmats_216nodes
list_dir = '.../lists/';  %path to the list with the participant IDs

addpath(bcn_dir)
addpath(data_dir)
addpath (list_dir)

ID = fileread(strcat(list_dir,'....txt'));
ID = strsplit(ID);
ID(end)=[]

data = {'ses-1_task-rest'};     

parcellation_code   = 'Schaefer200';
no_of_nodes         = 216;
wlet_scale          = 5; 
FisherZ             = 1;
method_of_thr       = 'rnk';
thr                 = (5:5:30)/100;
rep_modularity      = 1000;
tau                 = .75;
res_par             = 1;
K                   = 5; % Minimum network size (3 in Schaefer et al. 2018 17 networks)
file_suffix         = {'rnk5', 'rnk10', 'rnk15', 'rnk20', 'rnk25', 'rnk30'}; 
file_tag            = 'withprior_tau75';

Orig = load('origpart.csv') %entails which nodes belong to which network in the Schaefer atlas

%% Main loop

for k = 1:length(thr)

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
                [Y I]       = sort(abs(r_up), 'descend');
                index_thres = round(max(I)*thr(k));
                r_thres     = Y(index_thres);

                W(abs(W)<=abs(r_thres)) = 0;

            elseif strcmp(method_of_thr, 'abs') % Absolute thresholding
                W(abs(W)<thr(k)) = 0; 
            end
            A = (ones(no_of_nodes, no_of_nodes).*(W ~= 0)); % Binarize

            % Get the mean connectivity strength
            r_vec     = W(find(triu(W) & ~eye(no_of_nodes)));
            FC        = sum(r_vec)/length(r_vec);

            % Get the network metrics for the unweighted graph
            Eg       = efficiency_bin(A);
            Eloc     = efficiency_bin(A, 1);
            MeanEloc = sum(Eloc)/no_of_nodes;
            
            %Control random number generator
            rng(k*s*d)
            
            Qr = nan(1, rep_modularity); 
            Cr = nan(rep_modularity, no_of_nodes);
            
            for r = 1:rep_modularity
                [Cr(r, :), Qr(r)] = community_louvain(A, res_par, Orig);
            end

            [Qmax, i] = max(Qr);
            Q         = Qmax; 
            
            % Get the agreement matrix
            D = agreement(Cr');
    
            % Convert the elements of the agreement matrix into probabilities
            D = D./rep_modularity;
    
            % Apply consensus clustering across iterations, leaving the singleton
            % communities disconnected.
            Cc            = consensus_und(D, tau, rep_modularity);
            
            with_conn_norm = []; betw_conn_norm = []; ii=0; jj=0;
            for mm = 1:max(Cc)
                if length(find(Cc==mm))>K
                    ii = ii+1;
                    for mn = 1:max(Cc)
                        if length(find(Cc==mn))>K
                            if mn == mm
                                MM = W(find(Cc == mm), find(Cc == mn));
                                MM_up = MM(find(triu(MM) & ~eye(length(MM))));
                                with_conn_norm(ii) = sum(MM_up)/(size(MM, 1)*size(MM, 2)); % According to Garcia et al. (2018)
                            elseif mn>mm
                                jj = jj+1;
                                MN = W(find(Cc == mm), find(Cc == mn));
                                MN_vec = reshape(MN, size(MN, 1)*size(MN, 2), 1); 
                                betw_conn_norm(jj) = sum(MN_vec)/(length(find(Cc == mm))*length(find(Cc == mn)));
                            end
                        end
                    end
                end
            end
            
            %Create Co-classification matrix
            Cc1=Cc
            Cc2=Cc'
            Cc3=Cc2
            
            for i=1:(no_of_nodes-1) 
                Cc1=horzcat(Cc1,Cc);
                Cc3=vertcat(Cc3,Cc2);
            end
            
            CoCc= Cc1==Cc3

            fprintf('       Saving the network metrics ...\n')

            save(strcat(ID{s}, '_', data{d},'_scale', num2str(wlet_scale), '_', file_suffix{k}, '_', file_tag, '.mat'), ...
                'RawFC', 'Cc','CoCc','FC', 'Eg', 'Q', 'Eloc', 'MeanEloc', 'betw_conn_norm', 'with_conn_norm')

            fprintf('       DONE.\n')

        end
    end
end

cd(data_dir)

fprintf('   *** End of graph analysis. \n')
fprintf(['   ' datestr(now, 0) ' \n'])

%% End of the script