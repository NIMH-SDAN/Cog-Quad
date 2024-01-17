%% CO-CLASSIFICATION OF THE NETWORK MATRICES 
% Estimation of co-classification across all individuals   
% Produces the sample-specific networks
%                                                                         
% Requires functional connectivity maps.                                                                         
% Julia Linke July/2020                                               
%
%% Configuration
clear
clc

fprintf('   *** Averaging of network matrices per condition ...\n')

bcn_dir = ''; %path to the brain connectivity toolbox (download here: https://sites.google.com/site/bctnet/)
data_dir = '.../netmats'; %path to the output of the previous script b_Netmats_216nodes
list_dir = '.../lists/';  %path to the list with the participant IDs
out_dir = '.../CoClass/';

addpath(bcn_dir)
addpath(data_dir)
addpath(list_dir)

ID = fileread(strcat(list_dir,'....txt'));
ID = strsplit(ID);
ID(end)=[]

data = {'ses-1_task-rest'};

thr={'5','10','15','20','25','30'}
scale='5'
tau=.75
gtau='75'
tagtau='75'
rep_modularity = 1000;
sumArray=[]

cd(data_dir)
if ~exist(out_dir, 'dir')
    mkdir (out_dir)
end

%% Co-Classification
for d = 1:length(data)
    
    for t = 1:length(thr)

        for s = 1:length(ID)
            
            disp(['       Subject: ' ID{s} ', Threshold: ' thr{t} ', Data: ' data{d}])
            
            % Load the connectivity matrix
            cd (fullfile(data_dir,ID{s}))
            load(strcat(ID{s}, '_' ,data{d},'_scale',scale,'_rnk',thr{t},'_withprior_tau',tagtau,'.mat'));
            
            Cc1 = Cc;
            Cc2 = Cc';
            Cc3 = Cc2;
            
            for i=1:215 
                Cc1=horzcat(Cc1,Cc); 
                Cc3=vertcat(Cc3,Cc2);
            end
            
            CoCc = Cc1==Cc3
            
            if s==1
                sumArray = CoCc;
            else 
                sumArray = sumArray + CoCc;
            end
            clear Cc1 Cc2 Cc3 CoCc
        end
        
        meanArray   = sumArray / length(ID);
        meanArray   = meanArray - diag(diag(meanArray));
        sumArray   = sumArray - diag(diag(sumArray));
        
        
        %Consensus modularity based on co-classification matrix
        Ccm          = consensus_und(meanArray, tau, rep_modularity);
        Ccs          = consensus_und(sumArray, tau, rep_modularity);
       
        cd (out_dir)
        save(strcat('..._', data{d},'_scale', scale, '_rnk',thr{t},'_tau',tagtau,'_grouptau',gtau,'.mat'),'Ccm','Ccs','meanArray','sumArray','meanArray')
        
    end
    clear sumArray
    clear meanArray
end
    
cd(data_dir)

fprintf('   *** End of CoClassification. \n')
%% End of the script