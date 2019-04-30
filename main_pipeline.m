% Main pipeline for artifact removing
% Author: Albert J. Linan
%% 0 - It is needed to load eeglab and open a raw dataset
%% 1 - Creating the markers for the gradient artifact removing

% Inputs and outputs

    % Inputs
    
    % eeglab
    % alldata ~ ALLEEG struct
    % data ~ EEG struct

    % Time of repetition(ms)
    % TR = 1260; 

    % Threshold
    % This threshold should be selected loocking at the histogram but 200 works
    % TH = 200;

    % Figures
    % 1: Displays the histogram and the rest of the figures
    % 0: No 

    % s_channels 
    % It is a vector that contains the selected channels to calculate the 
    % average channel used for creating the markers

    % volumes
    % 1: Create volume/slice markers
    % 0: Create slices markers

% Outputs

    % new_alldata ~ ALLEEG struct
    % new_data ~ EEG struct
    % points ~ contains all the markers
tic;    
fprintf('\n--->>>--->EEG Artefact Correction starting...\n\n');
fprintf('\n--->>>---> 1 - Creating markers for the gradient artifact...\n\n');

TR = 1260;
TH = 200;
figures = 1;
s_channels = 1:31;
volumes = 1;

[ ~ , EEGMARK, puntos] = create_markers_v6(ALLEEG, EEG, TR, TH, figures, s_channels, volumes);

disp('Press a key to continue !') 
pause;

%% 2 - Removing the gradient artifact

fprintf('\n--->>>---> 2 - Removing the gradient artifact...\n\n');

% It is needed to cast the data to double
EEGMARK.data = double(EEGMARK.data);

% Load the functions from the following folder
addpath(genpath('fmrib1.21'))

% Apply artifact correction

lpf = [];               % Low pass filter cutoff (default: [ ]=70).
L = 1;                  % Interpolation folds (default: [ ]=10).
Win = 11;               %  Number of artifacts in avg. window (default: [ ]=30)
etype = 'Scan Start';   % Name of FMRI slice (slice timing) event.  Unlike
                        % fmrib_fastr.m, this is the name of the event.  fmrib_fastr.m 
                        % takes a vector with slice locations.
strig = 1;              % 1 for slice triggers (default)
                        % 0 for volume/section triggers
anc_chk = 0;            % 1 to perform ANC
                        % 0 No ANC
trig_correct = 0;       % 1 to correct for missing triggers;
                        % 0 to NOT correct.
Volumes = [];           % FMRI Volumes.  Needed if trig_correct=1; 
                        % otherwise use [ ];
Slices = [];            % FMRI Slices/Vol. usage like Volumes.
pre_frac = [];
exc_chan = 32;          
NPC = 'auto';           % Number of principal components to fit to residuals.
                        % 0 to skip OBS fitting and subtraction.
                        % 'auto' for automatic order selection (default)

[EEGGAR, command] = pop_fmrib_fastr(EEGMARK,lpf,L,Win,etype,...
    strig,anc_chk,trig_correct,Volumes,Slices,pre_frac,exc_chan,NPC); 

figure 
hold on
plot(EEGMARK.data(1,:))
plot(EEGGAR.data(1,:))

title(' Gradient Artifact Correction');
ylabel(' Amplitude(uV) ')
xlabel(' Points ')
legend('Original EEG', 'GA Removed EEG');
hold off

%% 3 - Downsampling the signal from 5000Hz to 250Hz

fprintf('\n--->>>---> 3 - Downsampling the signal from original Hz to 250Hz...\n\n');

freq = 250;
[EEGDWS] = pop_resample(EEGGAR, freq);

%% 4 - Creating the markers for the pulse artifact removing(qrs)

fprintf('\n--->>>---> 4 - Creating markers for the pulse artifact removing...\n\n');

ecgchan = 32;
[EEGQRS,~] = pop_fmrib_qrsdetect(EEGDWS,ecgchan,'qrs','no');

%% 5 - Applying ICA

fprintf('\n--->>>---> 5 - Applying ICA...\n\n');

% First is needed to change an option for saving the activations
pop_editoptions('option_computeica', 1 );
EEGICA = pop_runica( EEGQRS, 'icatype', 'runica', 'chanind', 1:31, 'extended', 1);

%% 6 - Applying EEG BCG Correction

fprintf('\n--->>>---> 6 - Applying BCG Correction...\n\n');

% Reference: https://github.com/rmabreu/BCG_Artefact_Correction

% change the following parameters accordingly
art_harm = 5;           % number of harmonics to include in the BCG artefact correction quantification
win_hz = 0.065;         % window length for which the BCG artefact correction will be assessed
filters = [ 0.5 45 ];   % low and high cutoff values of the band-pass filtering of EEG data
TR = TR/1000;           % repetition time of the fMRI acquisition [s]
harm_thr = 0.25;        % threshold for inclusion of harmonics {recommended = 0.25; for higher variability across subjects and channels, the threshold should be lower}
plt = 0;                % 1/0 = do/don't display the results
nb_plts = 4;            % number of subplots within each plot for displaying the selected IC time-courses


% load all necessary matrices (data, activations and mixing matrix from ICA)

dataset = EEGICA;
sph = EEGICA.icasphere;
wghts = EEGICA.icaweights;
activations = EEGICA.icaact;

mix_matrix = wghts * sph;

j = 1;

for i = 1:length(dataset.event)
    
   if strcmp(dataset.event(i).type,'qrs')
       
      Kp(j) = dataset.event(i).latency;  
      j = j + 1;
      
   end
      
end
    
dataset.Kp = unique(sort(round(Kp)));   
dataset.ecg = dataset.data(32,:);  
dataset.data = dataset.data(1:31,:);

% retrieve R peak annotations and compute appropriate limits for epoching
% EEG data using the R peaks as triggers (dependent on the subjects' heart
% rate)
R = dataset.Kp; 
m_R = (min(diff(R)) / dataset.srate) * 1000;
lim_inf_ms = -100;
if m_R < 4 * abs(lim_inf_ms)
    lim_sup_ms = lim_inf_ms + (4 * abs(lim_inf_ms));
else
    lim_sup_ms = lim_inf_ms + m_R;
end
limits_ecg = [ lim_inf_ms lim_sup_ms ];

% define the standard time-delay (210 ms) between R peak and BCG artefact
% occurrences (according to Allen et al., NeuroImage 1998)
delay_QRS = 0.21;

cont = 0;
while ~cont
 
    method = input('Choose the method Ballistocardiogram artifact correction: \n1 - PROJIC \n2 - PROJIC-OBS \n3 - PROJIC-AAS \n \n');
    cont = sum(method == [1 , 2 , 3 ]);
    
end

switch method
    
    case 1
        
        % PROJIC
        % Load the functions from the following folder
        addpath(genpath('BCG_Artefact_Correction-master'))
        % k_clusters = number of clusters. The algorithm can try with differents
        % numbers of cluster
        k_clusters = 2:5; 
        [ PowerFreq, PowerFreq_bkg ] = BCG_Correction_PROJIC(dataset, activations, ...
            mix_matrix, limits_ecg, filters, k_clusters, art_harm, harm_thr, TR, win_hz, ...
            plt, nb_plts);

        mttoc = floor(toc/60); sttoc = round(toc - mttoc*60);
        if mttoc < 60
            fprintf('--->>>--->EEG Artefact Correction Pipeline using PROJIC finished in %d min %d sec.\n', mttoc, sttoc);
        else
            httoc=floor(mttoc / 60); mttoc=round(mttoc - httoc*60);
            fprintf('--->>>--->EEG Artefact Correction Pipeline using PROJIC finished in %d hrs %d min %d sec.\n', httoc, mttoc, sttoc);
        end
        
    case 2
        
        %PROJIC-OBS
        % npc =  number of principal components
        k_clusters = 2:4; npc = 3:6; 
        [ PowerFreq_obs, PowerFreq_bkg_obs ] = BCG_Correction_PROJIC_OBS(dataset, activations, ...
            mix_matrix, limits_ecg, filters, npc, k_clusters, art_harm, harm_thr, ...
            TR, win_hz, delay_QRS, plt, nb_plts);

        mttoc = floor(toc/60); sttoc = round(toc - mttoc*60);
        if mttoc < 60
            fprintf('--->>>--->EEG Artefact Correction Pipeline using PROJIC-OBS finished in %d min %d sec.\n', mttoc, sttoc);
        else
            httoc=floor(mttoc / 60); mttoc=round(mttoc - httoc*60);
            fprintf('--->>>--->EEG Artefact Correction Pipeline using PROJIC-OBS finished in %d hrs %d min %d sec.\n', httoc, mttoc, sttoc);
        end

    case 3
        
        % PROJIC-AAS
        % n_win = numbers of averaging windows
        k_clusters = 2:5; n_win = 10:10:40; 
        [ PowerFreq_aas, PowerFreq_bkg_aas ] = BCG_Correction_PROJIC_AAS(dataset, activations, ...
            mix_matrix, limits_ecg, filters, n_win, k_clusters, art_harm, harm_thr, ...
            TR, win_hz, delay_QRS, plt, nb_plts);

        mttoc = floor(toc/60); sttoc = round(toc - mttoc*60);
        if mttoc < 60
            fprintf('--->>>--->EEG Artefact Correction Pipeline using PROJIC-AAS finished in %d min %d sec.\n', mttoc, sttoc);
        else
            httoc=floor(mttoc / 60); mttoc=round(mttoc - httoc*60);
            fprintf('--->>>--->EEG Artefact Correction Pipeline using PROJIC-AAS finished in %d hrs %d min %d sec.\n', httoc, mttoc, sttoc);
        end

end

% % Now the ICs related with the ECG are stored, the table has de following
% % structure =  [numclusters,n_components]
% ICclusters = cell2table(PowerFreq(3,:,:));
% [~,n_clusters] = size(ICclusters);
% ICclusters = ICclusters(:,2:n_clusters);
% ICclusters = table2array(ICclusters);


