function [ ] = run_uhdapter_audapter(outputdir,expt,bTestMode)
% outputdir: directory where data is saved
% expt: expt params set up in wrapper function
% h_fig: figure handles for display


if nargin < 1 || isempty(outputdir), outputdir = pwd, end
if nargin < 2, expt = []; end



%if nargin < 1, snum = []; end
%if nargin < 2, gender = []; end
if nargin < 3, bTestMode = 0; end
snum = [];
gender = [];


% set experiment-specific fields (or pass them in as 'expt')
expt.name = 'uhdapter';
if ~isfield(expt,'snum'), expt.snum = get_snum; end


if isnumeric(expt.snum)
    expt.snum = sprintf('s%02d',expt.snum);
elseif ~ischar(expt.snum)
    error('Subject ID must be a number or character string.')
end

expt.dataPath = get_exptSavePath('stress','acousticdata',expt.snum);


expt.gender = get_gender(gender);
stimtxtsize = 200;

% subjdir = fullfile(outputdir,'acousticdata',num2str(expt.snum))
if ~exist(expt.dataPath,'dir')
    mkdir(expt.dataPath)
end


trialdirname = 'temp_trials';
if bTestMode
    trialdir = fullfile(expt.dataPath,'practice',trialdirname);
else
    trialdir = fullfile(expt.dataPath,trialdirname);
end

if ~exist(trialdir,'dir')
    mkdir(trialdir)
end

%% set up stimuli
expt.conds = {'noisebase' 'baseline' 'ramp' 'hold' 'post' 'washout'}; % transfer task (to strings)?
expt.stress = {'trochee' 'iamb'};
expt.vowels = {'ey' 'eh' 'ah'};
expt.words = {'abate' 'adept' 'above' 'beta' 'meta'}; % data? bud?
%expt.noise = {'soft' 'mask'}


%% build up trial list

if ~isfield(expt,'startBlock')
    expt.startBlock = 1
end

expt.maxshift = 100

%% 1. alteration period
nNoisebase = 10*length(expt.words); %50
nBaseline = 22*length(expt.words); %110 to make numbers work out
nRamp = 4*length(expt.words); % 20
nHold = 50*length(expt.words); %250
nPost = 10*length(expt.words); %50
nWash = 4*length(expt.words);   %20

expt.nNoisebase = nNoisebase
expt.nBaseline=nBaseline
expt.nRamp=nRamp
expt.nHold=nHold
expt.nPost=nPost
expt.nWash=nWash

if bTestMode
    expt.maxshift = length(expt.words);
    nNoisebase = 2;
    nBaseline = 2;
    nRamp = length(expt.words);
    nHold = 2;
    nPost = 2;
    nWash = 2;
end

ramp = linspace(expt.maxshift./nRamp,expt.maxshift,nRamp); %ramped change in F1.
expt.shiftMags = [ones(1,nNoisebase) ones(1,nBaseline) ramp expt.maxshift.*ones(1,nHold) ones(1,nPost) ones(1,nWash)]; % build up list of formant shifts

expt.timing.stimdur = 1.75;                          % time stim is on screen, in seconds
expt.timing.interstimdur = .75;                   % minimum time between stims, in seconds
expt.timing.interstimjitter = .75;                % maximum extra time between stims (jitter)
expt.timing.visualfbdur = .25

expt.trial_dur = 1.75; % trial duration in s
expt.ntrials = nNoisebase + nBaseline + nRamp + nHold + nPost + nWash; %total trials
if bTestMode
    expt.ntrials_per_block = 5;
else
    expt.ntrials_per_block = 20;
end
expt.nblocks = expt.ntrials ./ expt.ntrials_per_block;
expt.allConds = [ones(1,nNoisebase) 2.*ones(1,nBaseline) 3.*ones(1,nRamp) 4.*ones(1,nHold) 5.*ones(1,nPost) 6.*ones(1,nWash)];
expt.listConds = expt.conds(expt.allConds); % vector with the condition indexed by the number in allConds

%% randomize by block

expt.listWords = []
expt.allWords = []
for rb = 1:expt.nblocks
    blockname = ['block_' num2str(rb)];
    block = repmat(1:length(expt.words),1,expt.ntrials_per_block/(length(expt.words)));
    allblockrp = block(randperm(expt.ntrials_per_block));
    expt.allWords = [expt.allWords allblockrp]
    listBlock = expt.words(allblockrp)
    expt.listWords = [expt.listWords listBlock]
end


expt.listVowels = txt2ipa(expt.listWords)

% many:1 word:vowel mapping
for t=1:expt.ntrials
    expt.allVowels(t) = find(strcmp(expt.listVowels{t},expt.vowels)); %get index of expt.vowels where a given vowel in listVowels matches
end

% list stress
expt.listStress = txt2mtr(expt.listWords);

% all stress
for t=1:expt.ntrials
    expt.allStress(t) = find(strcmp(expt.listStress{t},expt.stress)); %get index of expt.vowels where a given vowel in listVowels matches
end


% set noise

expt.listNoise = [2.*ones(1,(nNoisebase)) 3.*ones(1,(nBaseline+nRamp+nHold)) 2.*ones(1,nPost) 3.*ones(1,(nWash))] ;% 3 for speech + noise, 2 for just masking noise

% set missing expt fields to defaults
expt = set_exptDefaults(expt);

expt.inds = get_exptInds(expt,{'conds', 'words', 'vowels', 'stress'});

% save expt
if bTestMode
    save(fullfile(expt.dataPath,'practice','expt.mat'),'expt')
else
    save(fullfile(expt.dataPath,'expt.mat'), 'expt')
end

%% set up audapter
audioInterfaceName = 'Focusrite USB'; %SMNG default for Windows 10
Audapter('deviceName', audioInterfaceName);
Audapter('ost', '', 0);     % nullify online status tracking/
Audapter('pcf', '', 0);     % pert config files (use pert field instead)

% set audapter params
p = getAudapterDefaultParams(expt.gender); % get default params
% overwrite selected params with experiment-specific values:
p.bShift = 1;
p.bRatioShift = 0;
p.bMelShift = 1;



w = get_noiseSource(p);
Audapter('setParam', 'datapb', w, 1);
p.fb = 3;          % set feedback mode to 3: speech + noise; 2 = noise-masking feedback; 4 = speech-shaped. Noise has been measured at 77. 

p.fb3Gain = 0.014; % in current setup where headphones are plugged into socket 3 with knob at 3 o'clock
p.fb2Gain = 0.16; % 77



%% Initially initiate Audapter
AudapterIO('init',p);

%% run experiment
% setup figures
h_fig = setup_exptFigs;
get_figinds_audapter % names figs: stim = 1, ctrl = 2, dup = 3;
h_sub = get_subfigs_audapter(h_fig(ctrl))

% give instructions and wait for keypress
h_ready = draw_exptText(h_fig,.5,.5,expt.instruct.introtxt,expt.instruct.txtparams);
pause
delete_exptText(h_fig,h_ready)




% run trials

for iblock = expt.startBlock:expt.nblocks
    pause(1);
    if iblock == expt.startBlock
        startTrial = mod(expt.startTrial,expt.ntrials_per_block);
        if ~startTrial, startTrial = expt.ntrials_per_block; end
    else
        startTrial = 1;
    end
    
    
    for itrial = startTrial:expt.ntrials_per_block   % for each trial
        %         % pause if 'p' is pressed
        if get_pause_state(h_fig,'p')
            pause_trial(h_fig);
        end
        
        % set trial index
        trial_index = (iblock-1).* expt.ntrials_per_block + itrial;
        
        % plot trial number in experimenter view
        cla(h_sub(1))
        ctrltxt = sprintf('trial: %d/%d, cond: %s',trial_index,expt.ntrials,expt.listConds{trial_index});
        h_trialn = text(h_sub(1),0,0.5,ctrltxt,'Color','black', 'FontSize',30, 'HorizontalAlignment','center');
        
        % set new perturbation
        p.pertAmp = expt.shiftMags(trial_index) * ones(1, 257);
        p.pertPhi = expt.shiftAngles(trial_index) * ones(1, 257);
        Audapter('setParam','pertAmp',p.pertAmp)
        Audapter('setParam','pertPhi',p.pertPhi) % I seem to remember you can set multiple params in one line -- test this
        
        % set noise
        if p.fb ~= expt.listNoise(trial_index)
            p.fb = expt.listNoise(trial_index)
            Audapter('setParam', 'fb', p.fb)
            % initialize Audapter
            AudapterIO('init',p);
        end
        
        % run trial in Audapter
        Audapter('reset'); %reset Audapter
        Audapter('start'); %start trial
        
        % display stimulus
        txt2display = expt.listWords{trial_index};
        color2display = expt.colorvals{expt.allColors(trial_index)};
        h_text(1) = draw_exptText(h_fig,.5,.5,txt2display, 'Color',color2display, 'FontSize',stimtxtsize, 'HorizontalAlignment','center');
        pause(expt.timing.stimdur);
        
        % stop trial in Audapter
        Audapter('stop');
        
        % get data
        data = AudapterIO('getData');
 
        figure(h_fig(stim)) % change this later once it's fixed
        h_dur = rectangle%('Position',circ_pos,'Curvature',[1,1],'Facecolor', 'g');
        CloneFig(h_fig(stim),h_fig(dup))
        pause(expt.timing.interstimdur);
        
        % clear screen
        delete_exptText(h_fig,[h_text h_dur])
        clear h_text h_dur
        
        % plot shifted spectrogram
        figure(h_fig(ctrl))
        subplot(h_sub(2))
        show_spectrogram(data.signalIn, data.params.sr, 'noFig');
        tAxis = 0 : p.frameLen : p.frameLen * (size(data.fmts, 1) - 1);
        plot(tAxis/data.params.sr,data.fmts(:, 1 : 2), 'b');
        plot(tAxis/data.params.sr,data.sfmts(:, 1 : 2), 'g');
        
        % pause for viewing duration feedback
        pause(expt.timing.visualfbdur);
        
        % clear screen
        clear h_text h_dur
        
        % add intertrial interval
        pause(expt.timing.interstimdur);
        % add jitter
        pause(rand*expt.timing.interstimjitter)
        
        % save trial
        trialfile = fullfile(trialdir,[num2str(trial_index) '.mat']);
        save(trialfile,'data')
    end
    
    % end of block: display break text
    if iblock < expt.nblocks
        breaktext = sprintf('Time for a break!\n%d of %d trials done.\n\nPress the space bar to continue.',trial_index,expt.ntrials);
    else
        breaktext = sprintf('Thank you!\n\nPlease press the space bar to finish.');
    end
    h_break = draw_exptText(h_fig,.5,.5,breaktext,expt.instruct.txtparams);
    pause
    delete_exptText(h_fig,h_break)
    
end


%% write experiment data and metadata
alldata = struct;
for i = 1:expt.ntrials
    load(fullfile(trialdir,[num2str(i) '.mat']))
    names = fieldnames(data);
    for j = 1:length(names)
        alldata(i).(names{j}) = data.(names{j});
    end
end
clear data
data = alldata;

if bTestMode
    save(fullfile(expt.dataPath,'practice','data.mat'),'data')
else
    save(fullfile(expt.dataPath,'data.mat'), 'data')
end

% remove temp trial directory
fprintf('Removing temp directory... ')
rmdir(trialdir,'s');
fprintf('done.\n')
close(h_fig)