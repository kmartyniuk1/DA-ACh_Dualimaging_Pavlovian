# DA-ACh_Dualimaging_Pavlovian
% Close all figure, clear all variables
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% Define paths and data to analyze

path2data = '/Users/kellymartyniuk/Desktop/Kellendonk Lab/Projects/Dual imaging/Pavlovian/Analysis/Raw data/'; %Location of the dBASELINE_PERata
path2savefolder = '/Users/kellymartyniuk/Desktop/Kellendonk Lab/Projects/Dual imaging/Pavlovian/Analysis/Analyzed Data/'; %Path to save

%% Define the list of mice to analyzed

% Determine it by the folders included in the path2data
mice_list = dir(path2data);
for o = length(mice_list):-1:1
    if mice_list(o).isdir == 0
        mice_list(o) = [];
    else
        if strcmp(mice_list(o).name,'.') == 1 || strcmp(mice_list(o).name,'..') == 1
            mice_list(o) = [];
        end
    end
end
Nmice = length(mice_list);

%% Define sessions to exclude 
% For the moment we have to add them manually until we define an automatic 
% rejection method. A cell for each animal with the names of the sessions
Sessions2exclude{1} = []; % [] for adding all
Sessions2exclude{2} = [];
Sessions2exclude{3} = [];

%% Setup the variables for the data you want to extract

IdChannel = {'Grab','dLight'}; %Channels. Pos1=STREAM_STORE1 and Pos2=STREAM_STORE2
STREAM_STORE1 = {'A05A','B05B'}; % name of the 405 store
STREAM_STORE2 = {'A65A','B65B'}; % name of the 465 store
min2remove = 1; % Minutes to remove from the beginning of the recording
TRANGE = [-10 15]; % window size [start time relative to epoc onset, window duration]
BASELINE_PER = [-5 0]; % baseline period within our window
ARTIFACT = Inf; % optionally set an artifact rejection level
THRESHOLD.dLight.peak = 2; % Number of standar deviations to use for the thresholding
THRESHOLD.dLight.dip = 2; % Number of standar deviations to use for the thresholding
THRESHOLD.Grab.dip = 1; % Number of standar deviations to use for the thresholding
THRESHOLD.Grab.peak = 2; % Number of standar deviations to use for the thresholding
Max_lag = 2; % Max time lag in seconds to test for correlation
corr_type = 'Pearson'; % Pearson or Spearman
show_plot = 0; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
limits2plot.dLight = []; %If empty, automatically adjusted to the data. 
limits2plot.Grab = []; %If empty, automatically adjusted to the data.
reanalysis = 1; % If 1, the code runs for sessions already analyze. If 0, session excluded if analysed
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data
%% Loop for all the mice
for m = 1:Nmice
    % Define the path to the mouse and find the folders to analyze:
    path2mouse = [path2data,mice_list(m).name,'/'];
    sessions = dir(path2mouse);
    for o = length(sessions):-1:1
        if contains(sessions(o).name,mice_list(m).name) == 0
            sessions(o)=[];
        end
    end
    
    % Create a folder to save the data for this mouse and define the path
    if exist([path2savefolder,mice_list(m).name],'dir') == 0
        mkdir([path2savefolder,mice_list(m).name])
    end
    path2save_mouse = [path2savefolder,mice_list(m).name,'/'];
    
    behav_data2load = dir(path2mouse);
    for o = length(behav_data2load):-1:1
        if behav_data2load(o).isdir == 1
            behav_data2load(o) = [];
        end
    end

    %% Loop for all the sessions for the mouse
     for s = 1:length(sessions)
        IdChannel = {'Grab','dLight'};
        STREAM_STORE1 = {'A05A','B05B'}; % name of the 405 store
        STREAM_STORE2 = {'A65A','B65B'}; % name of the 465 store
        
        % Define the path to the session and create folder to save if needed:
        PATH2SESSION = [path2mouse,sessions(s).name];
        if exist([path2save_mouse,sessions(s).name],'dir') == 0
            mkdir([path2save_mouse,sessions(s).name])
        end
        PATH2SAVE = [path2save_mouse,sessions(s).name,'/'];
        
        % Check if results are already saved for this session
        done = exist([PATH2SAVE,'ITI_analysis.mat'],'file');
        if done == 0 || reanalysis == 1
            if exist([PATH2SAVE,'figures'],'dir') == 0
                mkdir([PATH2SAVE,'figures'])
            end
        end
            
            %% Read the data
            % Now read the specified data from our block into a Matlab structure.
            data = TDTbin2mat(PATH2SESSION, 'TYPE', {'epocs', 'scalars', 'streams'});
            
            % Get the trig times of the stimulus
            stim_times = [data.epocs.PC0_.onset data.epocs.PC0_.offset];
            %% Load the operant box information to align the time:
            % Identify the data from this session (only )
            idx = strfind(sessions(s).name,'_');
            session = sessions(s).name(idx(2)+1:end);
            data2load = [];
            for ii = 1:length(behav_data2load)
                if contains(behav_data2load(ii).name,session) == 1
                    data2load = ii;
                end
            end
            
           
            if ~isempty(data2load)
                load([path2mouse,behav_data2load(data2load).name],'Var','rawlist')
                
                for nm = 1:length(rawlist.name)
                    if contains(rawlist.name{nm},mice_list(m).name(5:end)) == 1
                        behav_data = nm;
                    end
                end
                
                Cplus_On_Rows = find(Var{behav_data}(:,2) == 0071);
                Cminus_On_Rows = find(Var{behav_data}(:,2) == 0061);
                Cplus_Off_Rows = find(Var{behav_data}(:,2) == 0072);
                Cminus_Off_Rows = find(Var{behav_data}(:,2) == 0062);
                t_cplus = [Var{behav_data}(Cplus_On_Rows,1) Var{behav_data}(Cplus_Off_Rows,1)];
                t_cminus = [Var{behav_data}(Cminus_On_Rows,1) Var{behav_data}(Cminus_Off_Rows,1)];
                t_all_stim = [t_cplus;t_cminus];
                [~,idx] = sort(t_all_stim(:,1));
                t_all_stim = t_all_stim(idx,:);
                
                % Obtain the linear function for the increasing mismatch
                % between the fiber photometry and the operant boxes systems
                x = t_all_stim(:,1);
                y = stim_times(:,1) - t_all_stim(:,1);
                P = polyfit(x,y,1);
                t_diff = P(1)*Var{behav_data}(:,1) + P(2);
                
                % Adjust the time events from the operant box to the corrected
                % time aligned to the fiber photometry
                Adjusted_Var = Var{behav_data};
                Adjusted_Var(:,1) = Adjusted_Var(:,1) + t_diff;
                
                % Get the new timing for the head entries and stims
                Head_entries_Rows = find(Adjusted_Var(:,2) == 1011);
                t_cplus = [Adjusted_Var(Cplus_On_Rows,1) Adjusted_Var(Cplus_Off_Rows,1)];
                t_cminus = [Adjusted_Var(Cminus_On_Rows,1) Adjusted_Var(Cminus_Off_Rows,1)];
                t_head_entry = Adjusted_Var(Head_entries_Rows,1);
                
                t_all_stim = [t_cplus;t_cminus];
                [~,idx] = sort(t_all_stim(:,1));
                t_all_stim = t_all_stim(idx,:);
                
                % Reject head entries happenning within half a second
                A = diff(t_head_entry);
                B = find(A < 0.5);
                t_head_entry(B+1) = [];
                  
                %Adjust the lengths of both channels
                minLength = min([length(data.streams.(STREAM_STORE1{1}).data) ...
                    length(data.streams.(STREAM_STORE1{2}).data) ...
                    length(data.streams.(STREAM_STORE2{1}).data)...
                    length(data.streams.(STREAM_STORE2{2}).data)]);
                
                for channel = 1:length(STREAM_STORE1)
                    if length(data.streams.(STREAM_STORE1{channel}).data) > minLength
                        data.streams.(STREAM_STORE1{channel}).data...
                            (minLength+1:length(data.streams.(STREAM_STORE1{channel}).data)) = [];
                        %(1:length(data.streams.(STREAM_STORE1{channel}).data)-minLength) = [];
                        
                    end
                    if length(data.streams.(STREAM_STORE2{channel}).data) > minLength
                        data.streams.(STREAM_STORE2{channel}).data...
                            (minLength+1:length(data.streams.(STREAM_STORE2{channel}).data)) = [];
                        %                     (1:length(data.streams.(STREAM_STORE2{channel}).data)-minLength) = [];
                    end
                end
                
                % Gets time vector for our data
                Fs = data.streams.(STREAM_STORE1{1}).fs;
                Max_idx = minLength;
                dt = 1/Fs;
                time = 0:dt:dt*(Max_idx-1);
                
                dummie = 1:length(time);
                dummie = dummie(time > 5*60);
                idx_5min = dummie(1);
                clear dummie
                
                time_dwn = downsample(time(idx_5min:end),10);
                N = length(time_dwn);
                
                clear time
                
                %%  %% Get the dFF for both sensors
                % Filter out the slow fluctuation
                ftype = 'high';
                n = 2;
                Wn = 0.05/((Fs/10)/2);
                [a,b] = butter(n,Wn,ftype);
                
                for channel = 1:length(STREAM_STORE1)
                    
                    % downsample 10x
                    data405 = data.streams.(STREAM_STORE1{channel}).data(idx_5min:end);
                    data405 = downsample(data405,10);
                    
                    data465 = data.streams.(STREAM_STORE2{channel}).data(idx_5min:end);
                    data465 = downsample(data465,10);
                    
                    %Baseline the 465 using the 405
                    bls = polyfit(data405, data465, 1);
                    Y_fit = bls(1) .* data405 + bls(2);
                    Y_dF = data465 - Y_fit;
                    
                    %dFF using 405 fit as baseline
                    dFF.(IdChannel{channel}).raw = 100*(Y_dF)./Y_fit;
                    dFF.(IdChannel{channel}).filt = filtfilt(a,b,double(dFF.(IdChannel{channel}).raw));
                end
                
                % Visualize the whole session:
                plotAllSession = 1; % If 0 don't plot
                for ooo = 1
                    if plotAllSession
                        figure
                        plot(time_dwn,dFF.dLight.filt,'b')
                        hold on
                        plot(time_dwn,dFF.Grab.filt,'r')
                        for o = 1:size(t_cplus,1)
                            if o == 1
                                plot([t_cplus(o,1) t_cplus(o,1)],[-4 12],'m')
                                plot([t_cminus(o,1) t_cminus(o,1)],[-4 12],'c')
                                plot([t_cplus(o,2) t_cplus(o,2)],[-4 12],'k')
                            else
                                plot([t_cplus(o,1) t_cplus(o,1)],[-4 12],'m','HandleVisibility','off')
                                plot([t_cplus(o,2) t_cplus(o,2)],[-4 12],'k','HandleVisibility','off')
                                plot([t_cminus(o,1) t_cminus(o,1)],[-4 12],'c','HandleVisibility','off')
                            end
                            plot([t_cminus(o,2) t_cminus(o,2)],[-4 12],'k','HandleVisibility','off')
                        end
                        plot(t_head_entry,ones(length(t_head_entry),1)*10,'*k')
                        xlim([time_dwn(1) time_dwn(end)])
                        xlabel('Time (s)')
                        ylim([-4 12])
                        ylabel('Signal size (A.U.)')
                        legend('dLight','Grab','CS+','CS-','Dipper','Head entries','Location','northoutside','NumColumns',6)
                    end
                end
                
                %% Trials analysis
                % Define time index for the trial period
                % Define time index for the trial period
                temp_t = time_dwn - time_dwn(1);
                dummie = 1:length(temp_t);
                dummie = dummie(temp_t >= abs(TRANGE(1)));
                idx_Init = dummie(1);
                dummie = 1:length(temp_t);
                dummie = dummie(temp_t >= abs(TRANGE(2)));
                idx_End = dummie(1);
                
                % Define time index for baseline correction
                dummie = 1:length(temp_t);
                dummie = dummie(temp_t >= abs(BASELINE_PER(1)));
                baseline_idx1 = dummie(1);
                dummie = 1:length(temp_t);
                dummie = dummie(temp_t >= abs(BASELINE_PER(2)));
                baseline_idx2 = dummie(1);
                
                n = idx_Init + idx_End; % Length of each trial
                t_trials = temp_t(1:n) - temp_t(idx_Init);
                
                
                Stim_data.CSplus.idx = ones(size(t_cplus,1),n)*nan;
                Stim_data.CSminus.idx = ones(size(t_cminus,1),n)*nan;
                Stim_data.CSplus.head_entries.idx = cell(size(t_cplus,1),1);
                Stim_data.CSminus.head_entries.idx = cell(size(t_cminus,1),1);
                Stim_data.CSplus.head_entries.times = cell(size(t_cplus,1),1);
                Stim_data.CSminus.head_entries.times = cell(size(t_cminus,1),1);
                for i = 1:length(IdChannel)
                    Stim_data.CSplus.dFF.(IdChannel{i}).raw = ones(size(t_cplus,1),n)*nan;
                    Stim_data.CSminus.dFF.(IdChannel{i}).raw = ones(size(t_cminus,1),n)*nan;
                    Stim_data.CSplus.dFF.(IdChannel{i}).baseline_corrected = ones(size(t_cplus,1),n)*nan;
                    Stim_data.CSminus.dFF.(IdChannel{i}).baseline_corrected = ones(size(t_cminus,1),n)*nan;
                end
                
                % CS+
                for o = 1:size(t_cplus,1)
                    ix = find(abs(time_dwn-t_cplus(o)) == min(abs(time_dwn-t_cplus(o))));
                    if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                        tmp = ix - (idx_Init-1):ix + idx_End;
                        Stim_data.CSplus.idx(o,:) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.CSplus.dFF.(IdChannel{i}).raw(o,:) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                        end
                    elseif ix+idx_End > length(time_dwn)
                        tmp = ix - (idx_Init-1):length(time_dwn);
                        Stim_data.CSplus.idx(o,1:length(tmp)) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.CSplus.dFF.(IdChannel{i}).raw(o,1:length(tmp)) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                        end
                    elseif ix-(idx_Init-1) < 1
                        tmp = 1:ix + idx_End;
                        Stim_data.CSplus.idx(o,(n-length(tmp)+1):length(tmp)) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.CSplus.dFF.(IdChannel{i}).raw(o,(n-length(tmp)+1):length(tmp)) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                        end
                    end
                    if ~isempty(data2load)
                        dummie = 1:length(t_head_entry);
                        Stim_data.CSplus.head_entries.idx{o} = dummie...
                            (t_head_entry >= time_dwn(tmp(1)) & t_head_entry <= time_dwn(tmp(end)));
                        Stim_data.CSplus.head_entries.times{o} = ...
                            t_head_entry(Stim_data.CSplus.head_entries.idx{o})-t_cplus(o);
                    end
                end
                for i = 1:length(IdChannel)
                    Stim_data.CSplus.dFF.(IdChannel{i}).baseline_corrected = ...
                        Stim_data.CSplus.dFF.(IdChannel{i}).raw - nanmean...
                        (Stim_data.CSplus.dFF.(IdChannel{i}).raw...
                        (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                end
                
                % CS-
                for o = 1:size(t_cminus,1)
                    ix = find(abs(time_dwn-t_cminus(o)) == min(abs(time_dwn-t_cminus(o))));
                    if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                        tmp = ix - (idx_Init-1):ix + idx_End;
                        Stim_data.CSminus.idx(o,:) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.CSminus.dFF.(IdChannel{i}).raw(o,:) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                        end
                    elseif ix+idx_End > length(time_dwn)
                        tmp = ix - (idx_Init-1):length(time_dwn);
                        Stim_data.CSminus.idx(o,1:length(tmp)) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.CSminus.dFF.(IdChannel{i}).raw(o,1:length(tmp)) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                        end
                    elseif ix-(idx_Init-1) < 1
                        tmp = 1:ix + idx_End;
                        Stim_data.CSminus.idx(o,(n-length(tmp)+1):length(tmp)) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.CSminus.dFF.(IdChannel{i}).raw(o,(n-length(tmp)+1):length(tmp)) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                        end
                    end
                    if ~isempty(data2load)
                        dummie = 1:length(t_head_entry);
                        Stim_data.CSminus.head_entries.idx{o} = dummie...
                            (t_head_entry >= time_dwn(tmp(1)) & t_head_entry <= time_dwn(tmp(end)));
                        Stim_data.CSminus.head_entries.times{o} = ...
                            t_head_entry(Stim_data.CSminus.head_entries.idx{o})-t_cminus(o);
                    end
                end
                for i = 1:length(IdChannel)
                    Stim_data.CSminus.dFF.(IdChannel{i}).baseline_corrected = ...
                        Stim_data.CSminus.dFF.(IdChannel{i}).raw - nanmean...
                        (Stim_data.CSminus.dFF.(IdChannel{i}).raw...
                        (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                end
                clear dummie temp_t
            %% Find mean stim data
            CSplus_Grab_mean= mean(Stim_data.CSplus.dFF.Grab.baseline_corrected, 'omitnan');
            CSplus_dLight_mean = mean (Stim_data.CSplus.dFF.dLight.baseline_corrected, 'omitnan');
                
            CSminus_Grab_mean= mean(Stim_data.CSminus.dFF.Grab.baseline_corrected, 'omitnan');
            CSminus_dLight_mean = mean (Stim_data.CSminus.dFF.dLight.baseline_corrected, 'omitnan');
                
                
            %% ITI Analysis
            % Identify the index of stimulation times to exclude them for the peak
            % selection
            if ~isempty(data2load)
                idxStim = unique(sort([reshape(Stim_data.CSplus.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.CSplus.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.CSplus.idx(:,idx_Init+1:end),2)),...
                reshape(Stim_data.CSminus.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.CSminus.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.CSminus.idx(:,idx_Init+1:end),2))]));
            else
                idxStim = unique(sort([reshape(Stim_data.CSplus.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.CSplus.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.CSplus.idx(:,idx_Init+1:end),2)),...
                    reshape(Stim_data.CSminus.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.CSminus.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.CSminus.idx(:,idx_Init+1:end),2))]));

            end
            
            idxStim_CSplus = Stim_data.CSplus.idx(:,idx_Init+1:end);
            idxStim_CSplus = reshape(idxStim_CSplus',1,size(idxStim_CSplus,1)*size(idxStim_CSplus,2));
            
            idxStim_CSminus = Stim_data.CSminus.idx(:,idx_Init+1:end);
            idxStim_CSminus = reshape(idxStim_CSminus',1,size(idxStim_CSminus,1)*size(idxStim_CSminus,2));
            
            
            dummie = 1:N;
            idx2include = setdiff(dummie,idxStim);
            clear dummie
            for i = 1:length(IdChannel)
                mean_dFF.(IdChannel{i}).raw = mean(double(dFF.(IdChannel{i}).filt(idx2include)));
                std_dFF.(IdChannel{i}).raw = std(double(dFF.(IdChannel{i}).filt(idx2include)));
            end
                
                %% Get overall correlation
                %%% Using raw data (filt version)
                signal1 = dFF.dLight.filt(idx2include);
                signal2 = dFF.Grab.filt(idx2include);
                save_name = 'Raw overall correlation with Grab lag';
                
                % Compute the correlation when lagging Grab agaist dLight
                [lag2plot,Ovrl_corr.raw] = ovrl_corr_calculation(signal1,signal2,...
                    time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
                clear signal1 signal2
                
                %% Using zscore data (filt version)
                signal1 = dFF.dLight.zscorefilt;
                signal2 = dFF.Grab.zscorefilt;
                save_name = 'Zscore overall correlation with Grab lag';
                
                Compute the correlation when lagging Grab agaist dLight
                [~,Ovrl_corr.zscore] = ovrl_corr_calculation(signal1,signal2,...
                    time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
                
                clear signal1 signal2
                close all
            %% Get correlation during the ITI
                signal1 = dFF.dLight.filt(idx2include);
                signal2 = dFF.Grab.filt(idx2include);
                save_name = 'Overall ITI correlation with Grab lag';
                
                [lag2plot,ITI_corr] = ovrl_corr_calculation(signal1,signal2,...
                    time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
                
                clear signal1 signal2
                close all
                %% Get overall correlation during trials only
                %%% Using raw data (filt version)
                
                %CSplus
                signal1 = dFF.dLight.filt(idxStim_CSplus);
                signal2 = dFF.Grab.filt(idxStim_CSplus);
                save_name = 'Raw ITI correlation with Grab lag';
                
                
                % Compute the correlation when lagging Grab agaist dLight
                [~,CSpluscorr.raw] = ovrl_corr_calculation(signal1,signal2,...
                    time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
                clear signal1 signal2
                
                %CSminus
                signal1 = dFF.dLight.filt(idxStim_CSminus);
                signal2 = dFF.Grab.filt(idxStim_CSminus);
                save_name = 'Raw ITI correlation with Grab lag';
                
                
                % Compute the correlation when lagging Grab agaist dLight
                [~,CSminuscorr.raw] = ovrl_corr_calculation(signal1,signal2,...
                    time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
                clear signal1 signal2
                         
%% %% Save
            if done == 0 || overwrite == 1
                Params.MICEId = mice_list(m).name;
                Params.SESSIONId = sessions(s).name;
                Params.THRESHOLD = THRESHOLD;
                
                save([PATH2SAVE,'Session_analysis.mat'],'Stim_data', 'stim_times','t_trials',...
                 'ITI_corr','CSminuscorr','CSpluscorr','lag2plot','CSplus_Grab_mean','CSplus_dLight_mean',...
                 'CSminus_Grab_mean','CSminus_dLight_mean')
                
            end
        end
     end
end

                
