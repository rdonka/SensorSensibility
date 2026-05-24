%% SENSOR SENSIBILITY SIGNAL PROCESSING
% This script contains code used to process raw fiber photometry recordings 
% to determine the acute effects of morphine administration across escalating 
% doses on DA activity assayed via in vivo fiber photometry. We compared 
% three sensors: VTA DA cell body GCaMP6f, NAcLS dLight1.3b, and NAcLS GRABDA2h.

% This script processes the raw photometry data, identifies and quantifies
% spontaneous transient events, extracts whole signal mean and MAD, and
% outputs processed data to csv files for statistical analysis in R.
% The signal processing was performed using PASTa (v1.1.0), an opensource 
% MATLAB toolbox for the processing of fiber photometry data and the 
% detection and quantification of spontaneous transient events. 
% See Donka et al 2025 for analysis details:

%% #1 Prepare the MATLAB environment
% Set up user path inputs
rootdirectory = getenv("USERPROFILE") + "\";
datapath = append(rootdirectory,'Box\JRoit Lab\Rachel\RD Data\Morphine Fiber Photometry\Morphine FP Dose Response\Standard Analysis\'); % Path for analysis files - this is where the data and keys are saved

analysisfolder = append(rootdirectory,'Box\JRoit Lab\Rachel\RD Analysis\Morphine Fiber Photometry\Acute Morphine Dose Response\Standard Analysis\'); % Folder to output analysis csv files to
subjectfigurepath = append(rootdirectory,'Box\JRoit Lab\Rachel\RD Analysis\Morphine Fiber Photometry\Acute Morphine Dose Response\Standard Analysis\Subject Figures\'); % Folder to output figures to
overallfigurepath = append(rootdirectory,'Box\JRoit Lab\Rachel\RD Analysis\Morphine Fiber Photometry\Acute Morphine Dose Response\Standard Analysis\Overall Figures\'); % Folder to output figures to
paperfigurepath = append(rootdirectory, 'Box\JRoit Lab\Publications\In Process\Morphine Fiber Photometry\Figures\Figure Panels\');

% Add data folders to MATLAB path
addpath(analysisfolder);
addpath(genpath(append(rootdirectory,'Box\JRoit Lab\Publications\In Process\Morphine Fiber Photometry\'))); % Path for GitHub repository

% Add PASTa repository folders to MATLAB path - only needed if PASTa is not installed as a MATLAB Toolbox 
% through the MATLAB Toolbox Add-On Explorer or the .mbtlx file in the PASTa repository.
addpath(genpath(append(rootdirectory,'Desktop\GitHub_MyRepositories\PASTa\'))); % Path for GitHub repository

% Load in experiment key names - Subject Key and File Key
subjectkeyname = 'Subject Key - Morphine FP Dose Response.csv'; % Name of csv file containing subject information; set to '' if not using a Subject Key
filekeyname = 'File Key - Morphine FP Dose Response.csv';

% Load colors
figurecolors = readtable('MorphineFPDoseResponse_ColorKey.csv', 'Decimal',',', 'Delimiter',',');

S1color = char(figurecolors.HexCode(find(strcmp(figurecolors.Category,'treatcolors_raw') & strcmp(figurecolors.ColorVariableName,'1'))));
S3color = char(figurecolors.HexCode(find(strcmp(figurecolors.Category,'treatcolors_raw') & strcmp(figurecolors.ColorVariableName,'3'))));
S5color = char(figurecolors.HexCode(find(strcmp(figurecolors.Category,'treatcolors_raw') & strcmp(figurecolors.ColorVariableName,'5'))));
S7color = char(figurecolors.HexCode(find(strcmp(figurecolors.Category,'treatcolors_raw') & strcmp(figurecolors.ColorVariableName,'7'))));

M2color = char(figurecolors.HexCode(find(strcmp(figurecolors.Category,'treatcolors_raw') & strcmp(figurecolors.ColorVariableName,'2'))));
M4color = char(figurecolors.HexCode(find(strcmp(figurecolors.Category,'treatcolors_raw') & strcmp(figurecolors.ColorVariableName,'4'))));
M6color = char(figurecolors.HexCode(find(strcmp(figurecolors.Category,'treatcolors_raw') & strcmp(figurecolors.ColorVariableName,'6'))));
M8color = char(figurecolors.HexCode(find(strcmp(figurecolors.Category,'treatcolors_raw') & strcmp(figurecolors.ColorVariableName,'8'))));

% Prep other variables
viruses = {'GCaMP6f', 'GRABDA2h', 'dLight1.3b'};
doses = {'2.50', '5.00', '7.50', '10.00'};
injtypes = {'S','M'};

skipexistingfigs = 0;

%% Read in file key and subject key CSV files and prepare data frame
% Treatments are coded by number in filekey.
    % 1: 2.5ml Saline   2: 2.5mg Morphine
    % 3: 5.0ml Saline   4: 5.0mg Morphine
    % 5: 7.5ml Saline   6: 7.5mg Morphine
    % 7: 10.0ml Saline  8: 10.0mg Morphine

% Load in csv files into tables - use the loadKeys function
[experimentkey_raw] = createExperimentKey(rootdirectory, subjectkeyname, filekeyname);
% Remove Session Exclude
sessionincludeidxs = find(~strcmp({experimentkey_raw.SessionExclude},'EXCLUDE'));
[experimentkey] = experimentkey_raw(sessionincludeidxs);


%% #3  Basic Protocol 2.5 - Extract data
% Extract raw data from blocks using the function 'extractTDTdata' to extract raw data blocks.
% Set up inputs
sigstreamnames = {'x65A', '465A', 'x465A'}; % All names of signal streams across files
baqstreamnames = {'x05A', '405A', 'x405A'}; % All names of background streams across files
rawfolderpaths = string({experimentkey.RawFolderPath})'; % Create string array of raw folder paths
extractedfolderpaths = string({experimentkey.ExtractedFolderPath})'; % Create string array of extracted folder paths

extractTDTdata(rawfolderpaths,extractedfolderpaths,sigstreamnames,baqstreamnames,'skipexisting',1,'trim',5); % extract data

%% #4  Basic Protocol 2.6 - Load fiber photometry data into structure
% Load previously extracted data blocks and tie to experiment key. 
% Each block is loaded as a row in the data structure.
[rawdata] = loadKeydata(experimentkey); % Load data based on the experiment key into the structure 'rawdata'

%% Fix recording issues
% add marks for 425 2.5ml saline -  folder '425-230214-112724'
rawdata(find(strcmp({rawdata.BlockFolder},'425-230214-112724'))).injt(1) = floor(977*rawdata(find(strcmp({rawdata.BlockFolder},'425-230214-112724'))).fs);
rawdata(find(strcmp({rawdata.BlockFolder},'425-230214-112724'))).injt(2) = rawdata(find(strcmp({rawdata.BlockFolder},'425-230214-112724'))).injt(1)+floor(40*rawdata(find(strcmp({rawdata.BlockFolder},'425-230214-112724'))).fs);

% remove extra k pulses for 424 saline 5mg - BlockFolder '424-230302-074946'
rawdata(find(strcmp({rawdata.BlockFolder},'424-230302-074946'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'424-230302-074946'))).injt(1:2);

% delete mistake k pulse from 426 morphine 2.5mg - BlockFolder '426-230217-075055'
rawdata(find(strcmp({rawdata.BlockFolder},'426-230217-075055'))).injt(1) = [];

% add second k pulse for 426 morphine 7.5mg - BlockFolder '426-230307-085729'
rawdata(find(strcmp({rawdata.BlockFolder},'426-230307-085729'))).injt(2) = rawdata(find(strcmp({rawdata.BlockFolder},'426-230307-085729'))).injt(1)+floor(60*rawdata(find(strcmp({rawdata.BlockFolder},'426-230307-085729'))).fs);

% remove extra k pulses for 427 saline 2.5mg - BlockFolder '427-230216-122419'
rawdata(find(strcmp({rawdata.BlockFolder},'427-230216-122419'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'427-230216-122419'))).injt(1:2);

% delete mistake k pulse from 427 5mg morphine -  BlockFolder '427-230224-094530'
%rawdata(find(strcmp({rawdata.BlockFolder},'427-230224-094530'))).injt(1) = [];

% delete mistake k pulses from 429 saline - BlockFolder '429-230216-161008'
rawdata(find(strcmp({rawdata.BlockFolder},'429-230216-161008'))).injt(1:2) = [];

% add second k pulse for 429 morphine 2.5mg - BlockFolder '429-230217-123011'
rawdata(find(strcmp({rawdata.BlockFolder},'429-230217-123011'))).injt(2) = rawdata(find(strcmp({rawdata.BlockFolder},'429-230217-123011'))).injt(1)+floor(60*rawdata(find(strcmp({rawdata.BlockFolder},'429-230217-123011'))).fs);

% remove extra k pulses for 429 saline 5mg - BlockFolder '429-230223-111159'
rawdata(find(strcmp({rawdata.BlockFolder},'429-230223-111159'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'429-230223-111159'))).injt(1:2);

% remove extra k pulses for 434 saline 5mg - BlockFolder '434-230411-085422'
rawdata(find(strcmp({rawdata.BlockFolder},'434-230411-085422'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'434-230411-085422'))).injt(1:2);

% remove extra k pulses for 434 morphine 10mg - BlockFolder '434-230426-095402'
rawdata(find(strcmp({rawdata.BlockFolder},'434-230426-095402'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'434-230426-095402'))).injt(1:2);

% remove extra k pulses for 436 morphine 10mg - BlockFolder '436-230512-121007'
rawdata(find(strcmp({rawdata.BlockFolder},'436-230512-121007'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'436-230512-121007'))).injt(1:2);

% remove extra k pulses for 439 saline 2.5mg - BlockFolder '439-230512-083537'
rawdata(find(strcmp({rawdata.BlockFolder},'439-230512-083537'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'439-230512-083537'))).injt(1:2);

% remove extra k pulses for 440 saline 2.5mg - BlockFolder '440-230420-075732'
rawdata(find(strcmp({rawdata.BlockFolder},'440-230420-075732'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'440-230420-075732'))).injt(1:2);

% remove extra k pulses for 440 morphine 10mg - BlockFolder '440-230512-101508'
rawdata(find(strcmp({rawdata.BlockFolder},'440-230512-101508'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'440-230512-101508'))).injt(1:2);

% remove extra k pulses for 443 morphine 10mg - BlockFolder '443-230603-074814'
rawdata(find(strcmp({rawdata.BlockFolder},'443-230603-074814'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'443-230603-074814'))).injt(1:2);

% remove extra k pulses for 444 morphine 10mg - BlockFolder '444-230601-094406'
rawdata(find(strcmp({rawdata.BlockFolder},'444-230601-094406'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'444-230601-094406'))).injt(1:2);

% remove extra k pulses for 445 saline 10mg - BlockFolder '445-230601-114749'
rawdata(find(strcmp({rawdata.BlockFolder},'445-230601-114749'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'445-230601-114749'))).injt(1:2);

% remove extra k pulses for 449 saline 10mg - BlockFolder '449-230628-075924'
rawdata(find(strcmp({rawdata.BlockFolder},'449-230628-075924'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'449-230628-075924'))).injt(1:2);

% remove extra k pulses for 449 morphine 10mg - BlockFolder '449-230629-080423'
rawdata(find(strcmp({rawdata.BlockFolder},'449-230629-080423'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'449-230629-080423'))).injt(1:2);

% add second k pulse for 451 saline 10mg - BlockFolder '451-230627-091334'
rawdata(find(strcmp({rawdata.BlockFolder},'451-230627-091334'))).injt(2) = rawdata(find(strcmp({rawdata.BlockFolder},'451-230627-091334'))).injt(1)+floor(60*rawdata(find(strcmp({rawdata.BlockFolder},'451-230627-091334'))).fs);

% remove extra k pulses for 451 morphine 10mg - BlockFolder '451-230628-093901'
rawdata(find(strcmp({rawdata.BlockFolder},'451-230628-093901'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'451-230628-093901'))).injt(1:2);

% add second k pulse for 452 morphine 5mg - BlockFolder '452-230616-124523'
rawdata(find(strcmp({rawdata.BlockFolder},'452-230616-124523'))).injt(2) = rawdata(find(strcmp({rawdata.BlockFolder},'452-230616-124523'))).injt(1)+floor(60*rawdata(find(strcmp({rawdata.BlockFolder},'452-230616-124523'))).fs);

% remove extra k pulses for 452 morphine 10mg - BlockFolder '452-230628-112929'
rawdata(find(strcmp({rawdata.BlockFolder},'452-230628-112929'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'452-230628-112929'))).injt(1:2);

% remove extra k pulses for 476 morphine 10mg - BlockFolder '476-240312-081237'
rawdata(find(strcmp({rawdata.BlockFolder},'476-240312-081237'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'476-240312-081237'))).injt(1:2);

% add marks for 477 10ml saline -  folder '477-240318-112359'
rawdata(find(strcmp({rawdata.BlockFolder},'477-240318-112359'))).injt(1) = floor(15.5*60*rawdata(find(strcmp({rawdata.BlockFolder},'477-240318-112359'))).fs);
rawdata(find(strcmp({rawdata.BlockFolder},'477-240318-112359'))).injt(2) = floor(16*60*rawdata(find(strcmp({rawdata.BlockFolder},'477-240318-112359'))).fs);

% remove mistaken 1st k pulse for 491 saline 2.5mg - BlockFolder '491-240731-100504'
rawdata(find(strcmp({rawdata.BlockFolder},'491-240731-100504'))).injt = rawdata(find(strcmp({rawdata.BlockFolder},'491-240731-100504'))).injt(2:3);

% add second k pulse for 491 saline 10mg - BlockFolder '491-240827-100855'
rawdata(find(strcmp({rawdata.BlockFolder},'491-240827-100855'))).injt(2) = rawdata(find(strcmp({rawdata.BlockFolder},'491-240827-100855'))).injt(1) + ceil(rawdata(find(strcmp({rawdata.BlockFolder},'491-240827-100855'))).fs*30);

%% #6  Basic Protocol 2.8 - Crop fiber photometry data streams
% Crop data: remove pre and post experimental session samples.
% NOTE: Cropping is the only part of the pipeline that will alter the loaded data fields. 
% To ensure you only complete this step once per analysis, it is reccomended to input structure 
% 'rawdata' to the function and output a new structure 'data'.
preinjectionlength = 12; % set up the length of the baseline period - used by the normalization function
postinjectionlength = 60; % set up the length of the post injection period - used by the normalization function

% Find session start, session end, and new injection locations
for eachfile = 1:length(rawdata)
    sessionstart = rawdata(eachfile).injt(1) - (floor(preinjectionlength*60*rawdata(eachfile).fs));
    sessionend = rawdata(eachfile).injt(2) + (floor(postinjectionlength*60*rawdata(eachfile).fs));

    rawdata(eachfile).sessionstart = sessionstart;
    if length(rawdata(eachfile).sig) >= sessionend
        rawdata(eachfile).sessionend = sessionend;
    else
        disp(append('WARNING: Session short for file ', num2str(eachfile)))
        rawdata(eachfile).sessionend = length(rawdata(eachfile).sig);
    end
end

% Crop data
cropstartfieldname = 'sessionstart'; % name of field with session start index
cropendfieldname = 'sessionend'; % name of field with session end index
streamfieldnames = {'sig', 'baq'}; % which streams to crop
epocsfieldnames = {'injt','sess'}; % which epocs to adjust to maintain relative position - OPTIONAL INPUT

[croppeddata] = cropFPdata(rawdata,cropstartfieldname,cropendfieldname, streamfieldnames,'epocsfieldnames', epocsfieldnames); % Output cropped data into new structure called data

clear rawdata

%% #7  Basic Protocol 3.1 - Subtract and filter fiber photometry data
% % Subset data by sensor
[croppeddata_GCaMP6f] = croppeddata(strcmp({croppeddata.Virus},'GCaMP6f'));
[croppeddata_dLight] = croppeddata(strcmp({croppeddata.Virus},'dLight1.3b'));
[croppeddata_GRABDA2h] = croppeddata(strcmp({croppeddata.Virus},'GRABDA2h'));

% Subtract and filter data
sigfieldname = 'sig';
baqfieldname = 'baq';
fsfieldname = 'fs';

[subtracteddata_GCaMP6f] = subtractFPdata(croppeddata_GCaMP6f,sigfieldname,baqfieldname,fsfieldname,'baqscalingfreqmin',10,'baqscalingfreqmax',80); % adds sigsub (subtracted stream) and sigfilt (subtracted and filtered stream) to data frame
[subtracteddata_dLight] = subtractFPdata(croppeddata_dLight,sigfieldname,baqfieldname,fsfieldname,'baqscalingfreqmin',10,'baqscalingfreqmax',80); % adds sigsub (subtracted stream) and sigfilt (subtracted and filtered stream) to data frame
[subtracteddata_GRABDA2h] = subtractFPdata(croppeddata_GRABDA2h,sigfieldname,baqfieldname,fsfieldname,'baqscalingfreqmin',8,'baqscalingfreqmax',50); % adds sigsub (subtracted stream) and sigfilt (subtracted and filtered stream) to data frame

subtracteddata = [subtracteddata_GCaMP6f; subtracteddata_dLight; subtracteddata_GRABDA2h];

clear croppeddata

clear croppeddata_GCaMP6f
clear croppeddata_dLight
clear croppeddata_GRABDA2h

%% #8  Basic Protocol 3.2 - Plot stream traces
% Plot whole session streams for each file. 
% Use plotTraces to plot all raw traces - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
% Manually save plots to allow for customization (addition of injection start/stop lines)
outputfiletype = '.png';
for eachfile = 1:length(subtracteddata)
    maintitle = append(num2str(subtracteddata(eachfile).SubjectID),' (',subtracteddata(eachfile).Virus,') - ',subtracteddata(eachfile).InjType,' ', num2str(subtracteddata(eachfile).Dose)); % Create title string for current plot
    plotfilepath = append(subjectfigurepath,'Raw Stream Traces\SessionTraces_',num2str(subtracteddata(eachfile).SubjectID),'_',subtracteddata(eachfile).Virus,'_',subtracteddata(eachfile).InjType,'_',num2str(subtracteddata(eachfile).Dose),outputfiletype);

    if isfile(plotfilepath) == 0 || skipexistingfigs == 1
        alltraces = plotTraces(subtracteddata,eachfile,maintitle); % Save plot into object for customization
        for eachtile = 1:5 % Add injection start/stop lines to each stream tile
            nexttile(eachtile)
            xline(subtracteddata(eachfile).injt(1),'--','Injection','Color','#C40300','FontSize',6)
            xline(subtracteddata(eachfile).injt(2),'--','Color','#C40300','FontSize',8)
        end
    
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]); % Manually save the figure
        exportgraphics(gcf,plotfilepath,'Resolution',300)
    else
        continue;
    end
end

%% Plot whole session FFT power plots for each file
% Use plotFFTs to plot all frequency magnitude plots - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
for eachfile = 1:length(subtracteddata)
    maintitle = append(num2str(subtracteddata(eachfile).SubjectID),' (',subtracteddata(eachfile).Virus,') - ',subtracteddata(eachfile).InjType,' ', num2str(subtracteddata(eachfile).Dose)); % Create title string for current plot
    plotfilepath = append(subjectfigurepath,'FFTs\SessionFFTPower_',num2str(subtracteddata(eachfile).SubjectID),'_',subtracteddata(eachfile).Virus,'_',subtracteddata(eachfile).InjType,'_',num2str(subtracteddata(eachfile).Dose),'.png');

    if isfile(plotfilepath) == 0 || skipexistingfigs == 1
        currFFTpower = plotFFTpower(subtracteddata,eachfile,maintitle,'fs');
    
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]);
        exportgraphics(gcf,plotfilepath,'Resolution',300)
    else
        continue;
    end
end


%% SUBTRACTED DATA CLEANING - NOTCH MAJOR ARTIFACT
% 436 M 7.5
idx_436_M_Dose3 = find([subtracteddata.SubjectID]==436 & [subtracteddata.Dose]==7.5 & strcmp({subtracteddata.InjType},'M'));
subtracteddata(idx_436_M_Dose3).artifactnotchstart(1) = 1380380;
subtracteddata(idx_436_M_Dose3).artifactnotchend(1) = 1381010;

% 443 S 7.5
idx_443_M_Dose3 = find([subtracteddata.SubjectID]==443 & [subtracteddata.Dose]==7.5 & strcmp({subtracteddata.InjType},'S'));
subtracteddata(idx_443_M_Dose3).artifactnotchstart(1) = 3978360;
subtracteddata(idx_443_M_Dose3).artifactnotchend(1) = 3978680;

% 476 S 2.5
idx_476_S_Dose1 = find([subtracteddata.SubjectID]==476 & [subtracteddata.Dose]==2.5 & strcmp({subtracteddata.InjType},'S'));
subtracteddata(idx_443_M_Dose3).artifactnotchstart(1) = 1150780;
subtracteddata(idx_443_M_Dose3).artifactnotchend(1) = 1151280;

subtracteddata(idx_443_M_Dose3).artifactnotchstart(2) = 1484970;
subtracteddata(idx_443_M_Dose3).artifactnotchend(2) = 1485400;

subtracteddata(idx_443_M_Dose3).artifactnotchstart(3) = 4153270;
subtracteddata(idx_443_M_Dose3).artifactnotchend(3) = 4153590;

% 479 M 2.5
idx_479_M_Dose1 = find([subtracteddata.SubjectID]==479 & [subtracteddata.Dose]==2.5 & strcmp({subtracteddata.InjType},'M'));
subtracteddata(idx_479_M_Dose1).artifactnotchstart(1) = 1025590;
subtracteddata(idx_479_M_Dose1).artifactnotchend(1) = 1025950;

subtracteddata(idx_479_M_Dose1).artifactnotchstart(2) = 1187000;
subtracteddata(idx_479_M_Dose1).artifactnotchend(2) = 1187410;

subtracteddata(idx_479_M_Dose1).artifactnotchstart(3) = 1569850;
subtracteddata(idx_479_M_Dose1).artifactnotchend(3) = 1570180;

subtracteddata(idx_479_M_Dose1).artifactnotchstart(4) = 3520920;
subtracteddata(idx_479_M_Dose1).artifactnotchend(4) = 3521220;

% 479 M 5
idx_479_M_Dose2 = find([subtracteddata.SubjectID]==479 & [subtracteddata.Dose]==5 & strcmp({subtracteddata.InjType},'M'));
subtracteddata(idx_479_M_Dose2).artifactnotchstart(1) = 2821940;
subtracteddata(idx_479_M_Dose2).artifactnotchend(1) = 2822320;

% 483 M 5
idx_483_M_Dose2 = find([subtracteddata.SubjectID]==483 & [subtracteddata.Dose]==5 & strcmp({subtracteddata.InjType},'M'));
subtracteddata(idx_483_M_Dose2).artifactnotchstart(1) = 3343990;
subtracteddata(idx_483_M_Dose2).artifactnotchend(1) = 3344380;

% 483 M 10
idx_483_M_Dose3 = find([subtracteddata.SubjectID]==483 & [subtracteddata.Dose]==10 & strcmp({subtracteddata.InjType},'M'));
subtracteddata(idx_483_M_Dose3).artifactnotchstart(1) = 1542570;
subtracteddata(idx_483_M_Dose3).artifactnotchend(1) = 1543020;

% 483 S 5
idx_483_S_Dose2 = find([subtracteddata.SubjectID]==483 & [subtracteddata.Dose]==5 & strcmp({subtracteddata.InjType},'S'));
subtracteddata(idx_483_S_Dose2).artifactnotchstart(1) = 3268930;
subtracteddata(idx_483_S_Dose2).artifactnotchend(1) = 3269450;

%% NOTCH ARTIFACTS
streamfieldnames = {'sig','baq','baqscaled','sigsub','sigfilt'};
[data] = notchStreamArtifacts(subtracteddata,streamfieldnames,'artifactnotchstart','artifactnotchend');

%% #8  Basic Protocol 3.2 - Plot stream traces
% Plot whole session streams for each file. 
% Use plotTraces to plot all raw traces - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
% Manually save plots to allow for customization (addition of injection start/stop lines)
outputfiletype = '.png';
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).SubjectID),' (',data(eachfile).Virus,') - ',data(eachfile).InjType,' ', num2str(data(eachfile).Dose)); % Create title string for current plot
    plotfilepath = append(subjectfigurepath,'Raw Stream Traces\SessionTraces_SubtractedPostArtifactRemoval_',num2str(data(eachfile).SubjectID),'_',data(eachfile).Virus,'_',data(eachfile).InjType,'_',num2str(data(eachfile).Dose),outputfiletype);

    if isfile(plotfilepath) == 0 || skipexistingfigs == 1
        alltraces = plotTraces(data,eachfile,maintitle); % Save plot into object for customization
        for eachtile = 1:5 % Add injection start/stop lines to each stream tile
            nexttile(eachtile)
            xline(data(eachfile).injt(1),'--','Injection','Color','#C40300','FontSize',6)
            xline(data(eachfile).injt(2),'--','Color','#C40300','FontSize',8)
        end
    
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]); % Manually save the figure
        exportgraphics(gcf,plotfilepath,'Resolution',300)
    else
        continue;
    end
end

%% Trim out injection
streamfieldnames = {'sig','baq','baqscaled','sigsub','sigfilt'};

for eachfile = 1:length(data)
    for eachstream = 1:length(streamfieldnames)
        streamname = streamfieldnames{eachstream};
        disp(['Removing Injection: File ', num2str(eachfile)])
        allindices = (1:length(data(eachfile).(streamname)));
        includeindices = (allindices < data(eachfile).injt(1) | allindices > data(eachfile).injt(2));
        data(eachfile).(streamname) = data(eachfile).(streamname)(includeindices);
    end
end


%% Normalize subtracted and filtered data stream
% Normalize to session baseline mean and standard deviation
for eachfile = 1:length(data)
    data(eachfile).BLstart = 1;
    data(eachfile).BLend =  data(eachfile).injt(1)-1;
end

[data] = normBaseline(data,'sigfilt','BLstart','BLend');


%% Plot normalized streams
% Use plotNormTraces to plot all raw traces - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
streams = {'sigfiltz_normbaseline'};
streamtitles = {'Normalized to Baseline'};
outputfiletype = '.png';

for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).SubjectID),' (',data(eachfile).Virus,') - ',data(eachfile).InjType,' ', num2str(data(eachfile).Dose)); % Create title string for current plot    plotfilepath = append(figurefolder,'SessionBaselineNormZ_',num2str(data(eachfile).SubjectID),'_',data(eachfile).Phase,'_',num2str(data(eachfile).Session));
    plotfilepath = append(subjectfigurepath,'Normalized Stream Traces\SessionNormalizedTraces_',num2str(data(eachfile).SubjectID),'_',data(eachfile).Virus,'_',data(eachfile).InjType,'_',num2str(data(eachfile).Dose),outputfiletype);
    
    if isfile(plotfilepath) == 0 || skipexistingfigs == 1
        normtraces = plotNormTraces(data,eachfile,streams,'fs',maintitle,streamtitles);
    
        xline(data(eachfile).injt(1),'--','Injection','Color','#C40300','FontSize',8) 
    
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 10, 3]);
        exportgraphics(gcf,plotfilepath,'Resolution',300)
    else
        continue;
    end
end

%% #13 Basic Protocol 4.2 - Find and quantify transient events
for eachfile = 1:length(data)
    % Add thresholds for each session baseline
    data(eachfile).SDthreshold = 2.6;  
end

% Create list of variables to add to the new data structure with the output transient events
addvariablesfieldnames = [fieldnames(experimentkey); {'params'}]; % This makes a list of all fieldnames in the experimentkey and adds the 'params' field

% Find transients based on pre-peak baseline window mean - reccomended as the first pass choice for transient analysis
[transientdata_sigfiltznormBL] = findTransients(data,addvariablesfieldnames,'sigfiltz_normbaseline','SDthreshold','fs','AUCwindowms',1500);
[transientdata_sigfiltznormBL_AUCwindow1000] = findTransients(data,addvariablesfieldnames,'sigfiltz_normbaseline','SDthreshold','fs','AUCwindowms',1000);
[transientdata_sigfiltznormBL_AUCwindow2000] = findTransients(data,addvariablesfieldnames,'sigfiltz_normbaseline','SDthreshold','fs','AUCwindowms',2000);
[transientdata_sigfiltznormBL_AUCwindow3000] = findTransients(data,addvariablesfieldnames,'sigfiltz_normbaseline','SDthreshold','fs','AUCwindowms',3000);

% Bin transients with 3 minute bins
[transientdata_sigfiltznormBL] = binTransients(transientdata_sigfiltznormBL,'binlengthmins',3);
[transientdata_sigfiltznormBL_AUCwindow1000] = binTransients(transientdata_sigfiltznormBL_AUCwindow1000,'binlengthmins',3);
[transientdata_sigfiltznormBL_AUCwindow2000] = binTransients(transientdata_sigfiltznormBL_AUCwindow2000,'binlengthmins',3);
[transientdata_sigfiltznormBL_AUCwindow3000] = binTransients(transientdata_sigfiltznormBL_AUCwindow3000,'binlengthmins',3);

% Add Phase
for eachfile = 1:length(transientdata_sigfiltznormBL)
    currfiletransientquant = transientdata_sigfiltznormBL(eachfile).transientquantification;
    preinjidxs = currfiletransientquant.Bin_3min <= 4;
    postinjidxs = currfiletransientquant.Bin_3min > 4;

    currfiletransientquant.Phase = strings(height(currfiletransientquant),1);
    currfiletransientquant.Phase(preinjidxs) = "PREINJ";
    currfiletransientquant.Phase(postinjidxs) = "POSTINJ";

    transientdata_sigfiltznormBL(eachfile).transientquantification = currfiletransientquant;
end


%% #15 Basic Protocol 4.4 - Export transient events
% Saves all individual transient events to one table and exports the table to a csv file
addvariables = {'SubjectID','TreatNum','SessionExclude','Dose','InjType','Weight','InjVol','Sex','FiberPlacement','FiberSide','Virus','Rig','Power'};
alltransients_sigfiltznormBL = exportTransients(transientdata_sigfiltznormBL,'transientquantification',analysisfolder,addvariables,'exportfilename','transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold.csv');
alltransients_sigfiltznormBL_AUCwindow1000 = exportTransients(transientdata_sigfiltznormBL_AUCwindow1000,'transientquantification',analysisfolder,addvariables,'exportfilename','transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold_AUCwindow1000ms.csv');
alltransients_sigfiltznormBL_AUCwindow2000 = exportTransients(transientdata_sigfiltznormBL_AUCwindow2000,'transientquantification',analysisfolder,addvariables,'exportfilename','transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold_AUCwindow2000ms.csv');
alltransients_sigfiltznormBL_AUCwindow3000 = exportTransients(transientdata_sigfiltznormBL_AUCwindow3000,'transientquantification',analysisfolder,addvariables,'exportfilename','transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold_AUCwindow3000ms.csv');


%% Subset saline first session only
alltransients_sigfiltznormBL_S1idxs = find(strcmp(alltransients_sigfiltznormBL.InjType,'S')&[alltransients_sigfiltznormBL.Dose{:}]'==2.5);
alltransients_sigfiltznormBL_S1 = alltransients_sigfiltznormBL(alltransients_sigfiltznormBL_S1idxs,:);
writetable(alltransients_sigfiltznormBL_S1,append(analysisfolder,'transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold_S1Only.csv'))


%% ANALYZE TRANSIENT SHAPE
%% Plot color grade maps by bin for each file
% Use plotTransientColorGradeBins to plot color grade maps of all transient
% events by bin
transientstreamfieldname = 'transientstreamdatacentered';
binfieldname = 'Bin_3min';
outputfiletype = 'png';

for eachfile = 1:length(transientdata_sigfiltznormBL)
    maintitle = append(num2str(transientdata_sigfiltznormBL(eachfile).SubjectID),' (',transientdata_sigfiltznormBL(eachfile).Virus,') - ',transientdata_sigfiltznormBL(eachfile).InjType,' ',...
        num2str(transientdata_sigfiltznormBL(eachfile).Dose)); % Create title string for current plot
    plotfilepath = append(subjectfigurepath,'Transient Color Grade Maps\TransientBinColorGradeMap_LocalScaling_',num2str(transientdata_sigfiltznormBL(eachfile).SubjectID),'_',transientdata_sigfiltznormBL(eachfile).Virus,'_',...
        transientdata_sigfiltznormBL(eachfile).InjType,'_',num2str(transientdata_sigfiltznormBL(eachfile).Dose));

    currfiletransientcolormapbins = plotTransientColorGradeBins(transientdata_sigfiltznormBL,eachfile,transientstreamfieldname,binfieldname,maintitle,'climtype','local');

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 10, 8]);
    exportgraphics(gcf,append(plotfilepath,'.',outputfiletype),'Resolution',300)
end

%% Plot overall color grade maps by phase for each file
% Add Phase to transient table
for eachfile = 1:length(transientdata_sigfiltznormBL)
    preinjidxs = [transientdata_sigfiltznormBL(eachfile).transientquantification.Bin_3min] <= 4;
    postinjidxs = [transientdata_sigfiltznormBL(eachfile).transientquantification.Bin_3min] > 4;
    
    transientquant = transientdata_sigfiltznormBL(eachfile).transientquantification;

    transientquant.Phase = strings(height(transientquant),1);
    transientquant.Phase(preinjidxs) = 'PREINJ';
    transientquant.Phase(postinjidxs) = 'POSTINJ';
    
    transientdata_sigfiltznormBL(eachfile).transientquantification = transientquant;
end

% Use plotTransientColorGradePhases to plot color grade maps of all transient
% events by bin
transientstreamfieldname = 'transientstreamdatacentered';
phasefieldname = 'Phase';
outputfiletype = 'png';

for eachfile = 1:length(transientdata_sigfiltznormBL)
    maintitle = append(num2str(transientdata_sigfiltznormBL(eachfile).SubjectID),' (',transientdata_sigfiltznormBL(eachfile).Virus,') - ',transientdata_sigfiltznormBL(eachfile).InjType,' ',...
        num2str(transientdata_sigfiltznormBL(eachfile).Dose)); % Create title string for current plot
    plotfilepath = append(subjectfigurepath,'Transient Color Grade Maps\TransientPhaseColorGradeMap_LocalScaling_',num2str(transientdata_sigfiltznormBL(eachfile).SubjectID),'_',transientdata_sigfiltznormBL(eachfile).Virus,'_',...
        transientdata_sigfiltznormBL(eachfile).InjType,'_',num2str(transientdata_sigfiltznormBL(eachfile).Dose));

    currfiletransientcolormapphases = plotTransientColorGradePhases(transientdata_sigfiltznormBL,eachfile,transientstreamfieldname,phasefieldname,maintitle,'climtype','local');

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]);
    exportgraphics(gcf,append(plotfilepath,'.',outputfiletype),'Resolution',300)
end


%% Plot overall color grade maps by sensor and treatment
customColorScale = createCustomColorScale();
addvariablesfieldnames = {'SubjectID','TreatNum','SessionExclude','Dose','InjType','Weight','InjVol','Sex','FiberPlacement','FiberSide','Virus','Rig','Power'};
transientstreamfieldnames = {'transientstreamdata', 'transientstreamdatacentered'};

[alltransientstruct_sigfiltnormzBL] = exportTransientStreams(transientdata_sigfiltznormBL,transientstreamfieldnames,addvariablesfieldnames);

% Add Group to transient table
for eachevent = 1:height(alltransientstruct_sigfiltnormzBL)
    currdose = num2str(alltransientstruct_sigfiltnormzBL(eachevent).Dose);
    currinjtype = alltransientstruct_sigfiltnormzBL(eachevent).InjType;
    currphase = alltransientstruct_sigfiltnormzBL(eachevent).Phase;

    alltransientstruct_sigfiltnormzBL(eachevent).Group = append(currdose,'_',currinjtype,'_',currphase);
end

% Plot transient color maps by virus and treatment - uncentered - global clim scaling
outputfiletype = 'png';

for eachvirus = 1:length(viruses)
    currvirus = viruses{eachvirus};

    for eachdose = 1:length(doses)
        currdose = str2num(doses{eachdose});
        currfilepath = append(overallfigurepath,'Transient Color Grade Maps\TransientColorMap_AllEvents_sigfiltznormBL_global_',currvirus,'_',num2str(currdose));
        curralltransientstruct = alltransientstruct_sigfiltnormzBL(strcmp({alltransientstruct_sigfiltnormzBL.Virus},currvirus) & strcmp([alltransientstruct_sigfiltnormzBL.Phase],'POSTINJ') ...
            & [alltransientstruct_sigfiltnormzBL.Dose] == currdose);
        
        alltransientcolormaps = plotTransientColorGradeGroups(curralltransientstruct,'InjType','transientstreamdata',customColorScale,currvirus,'climtype','global');
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]);

        exportgraphics(gcf,append(currfilepath,'.',outputfiletype),'Resolution',300)
    end
end

% Plot transient color maps by virus and treatment - uncentered - local clim scaling
outputfiletype = 'png';

for eachvirus = 1:length(viruses)
    currvirus = viruses{eachvirus};

    for eachdose = 1:length(doses)
        currdose = str2num(doses{eachdose});
        currfilepath = append(overallfigurepath,'Transient Color Grade Maps\TransientColorMap_AllEvents_sigfiltznormBL_local_',currvirus,'_',num2str(currdose));
        curralltransientstruct = alltransientstruct_sigfiltnormzBL(strcmp({alltransientstruct_sigfiltnormzBL.Virus},currvirus) & strcmp([alltransientstruct_sigfiltnormzBL.Phase],'POSTINJ') ...
            & [alltransientstruct_sigfiltnormzBL.Dose] == currdose);
        
        alltransientcolormaps = plotTransientColorGradeGroups(curralltransientstruct,'InjType','transientstreamdata',customColorScale,currvirus,'climtype','local');
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]);

        exportgraphics(gcf,append(currfilepath,'.',outputfiletype),'Resolution',300)
    end
end


% Plot transient color maps by virus and treatment - centered - global clim scaling
outputfiletype = 'png';

for eachvirus = 1:length(viruses)
    currvirus = viruses{eachvirus};

    for eachdose = 1:length(doses)
        currdose = str2num(doses{eachdose});
        currfilepath = append(overallfigurepath,'Transient Color Grade Maps\TransientColorMap_AllEvents_sigfiltznormBL_centered_global_',currvirus,'_',num2str(currdose));
        curralltransientstruct = alltransientstruct_sigfiltnormzBL(strcmp({alltransientstruct_sigfiltnormzBL.Virus},currvirus) & strcmp([alltransientstruct_sigfiltnormzBL.Phase],'POSTINJ') ...
            & [alltransientstruct_sigfiltnormzBL.Dose] == currdose);
        
        alltransientcolormaps = plotTransientColorGradeGroups(curralltransientstruct,'InjType','transientstreamdatacentered',customColorScale,currvirus,'climtype','global');
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]);

        exportgraphics(gcf,append(currfilepath,'.',outputfiletype),'Resolution',300)
    end
end

% Plot transient color maps by virus and treatment - centered - local clim scaling
outputfiletype = 'png';

for eachvirus = 1:length(viruses)
    currvirus = viruses{eachvirus};

    for eachdose = 1:length(doses)
        currdose = str2num(doses{eachdose});
        currfilepath = append(overallfigurepath,'Transient Color Grade Maps\TransientColorMap_AllEvents_sigfiltznormBL_centered_local_',currvirus,'_',num2str(currdose));
        curralltransientstruct = alltransientstruct_sigfiltnormzBL(strcmp({alltransientstruct_sigfiltnormzBL.Virus},currvirus) & strcmp([alltransientstruct_sigfiltnormzBL.Phase],'POSTINJ') ...
            & [alltransientstruct_sigfiltnormzBL.Dose] == currdose);
        
        alltransientcolormaps = plotTransientColorGradeGroups(curralltransientstruct,'InjType','transientstreamdatacentered',customColorScale,currvirus,'climtype','local');
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]);

        exportgraphics(gcf,append(currfilepath,'.',outputfiletype),'Resolution',300)
    end
end


% Plot transient color maps by InjType and Dose - centered
for eachinjtype = 1:length(injtypes)
    currinjtype = injtypes{eachinjtype};
    for eachdose = 1:length(doses)
        currdose = str2num(doses{eachdose});
        currfilepath = append(overallfigurepath,'Transient Color Grade Maps\TransientColorMap_AllEvents_sigfiltznormBL_centered_',currinjtype,'_',num2str(currdose));
        currmaintitle = append('Dose: ',num2str(currdose),' - ',currinjtype);
        curralltransientstruct = alltransientstruct_sigfiltnormzBL(strcmp({alltransientstruct_sigfiltnormzBL.InjType},currinjtype) & strcmp([alltransientstruct_sigfiltnormzBL.Phase],'POSTINJ') ...
            & [alltransientstruct_sigfiltnormzBL.Dose] == currdose);
        
        alltransientcolormaps = plotTransientColorGradeGroups(curralltransientstruct,'Virus','transientstreamdatacentered',currmaintitle,'climtype','local','yscaletype','manual','yscalemin',-6,'yscalemax',12);
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 9, 3]);

        exportgraphics(gcf,append(currfilepath,'.',outputfiletype),'Resolution',300)
    end
end

% Plot transient color maps by InjType and Dose - centered
for eachinjtype = 1:length(injtypes)
    currinjtype = injtypes{eachinjtype};
    for eachdose = 1:length(doses)
        currdose = str2num(doses{eachdose});
        currfilepath = append(overallfigurepath,'Transient Color Grade Maps\TransientColorMap_AllEvents_sigfiltznormBL_centered_',currinjtype,'_',num2str(currdose));
        currmaintitle = append('Dose: ',num2str(currdose),' - ',currinjtype);
        curralltransientstruct = alltransientstruct_sigfiltnormzBL(strcmp({alltransientstruct_sigfiltnormzBL.InjType},currinjtype) & strcmp([alltransientstruct_sigfiltnormzBL.Phase],'POSTINJ') ...
            & [alltransientstruct_sigfiltnormzBL.Dose] == currdose);
        
        alltransientcolormaps = plotTransientColorGradeGroups(curralltransientstruct,'Virus','transientstreamdatacentered',currmaintitle,'climtype','local','yscaletype','manual','yscalemin',-6,'yscalemax',12);
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 9, 3]);

        exportgraphics(gcf,append(currfilepath,'.',outputfiletype),'Resolution',300)
    end
end


%% #16 Basic Protocol 4.5 - Plot whole session transients
% Use plotTransients to plot the whole session trace with detected transients for each file.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).SubjectID),' (',data(eachfile).Virus,') - ',data(eachfile).InjType,' ', num2str(data(eachfile).Dose)); % Create title string for current plot
    plotfilepath = append(subjectfigurepath,'Transients/SessionTransients_',num2str(data(eachfile).SubjectID),'_',data(eachfile).Virus,'_',data(eachfile).InjType,'_',num2str(data(eachfile).Dose));

    if isfile(plotfilepath) == 0 || skipexistingfigs == 1
        sessiontransients = plotTransients(data,eachfile,'sigfiltz_normbaseline','fs',transientdata_sigfiltznormBL,maintitle,'saveoutput',1,'outputfiletype','png','plotfilepath',plotfilepath);
    else
        continue;
    end
end

%% #17 Basic Protocol 4.6 - Plot binned transients
% Use plotTransientBins to plot session bin traces with detected transients for each file.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).SubjectID),' (',data(eachfile).Virus,') - ',data(eachfile).InjType,' ', num2str(data(eachfile).Dose)); % Create title string for current plot
    plotfilepath = append(subjectfigurepath,'Transients/SessionTransientsBins_',num2str(data(eachfile).SubjectID),'_',data(eachfile).Virus,'_',data(eachfile).InjType,'_', num2str(data(eachfile).Dose));
    if isfile(plotfilepath) == 0 || skipexistingfigs == 1
        allbins = plotTransientBins(data,eachfile,'sigfiltz_normbaseline','fs',transientdata_sigfiltznormBL,'Bin_3min',maintitle,'saveoutput',1,'outputfiletype','png','plotfilepath',plotfilepath);
    else
        continue;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PAPER FIGURES

%% Figure 2 - VTA GCaMP6f Example Morphine Traces
plotexampleSidx = 31; % 447 S 10
plotexampleMidx = 32; % 434 M 10

% Prep plot variables
currxlength = length(data(plotexampleMidx).sigfiltz_normbaseline);    
currxticklabels = 0:6:72;
currxticks = floor(currxticklabels.*60.*data(plotexampleMidx).fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz = max(data(plotexampleMidx).sigfiltz_normbaseline);
yminsigfiltz = min(data(plotexampleMidx).sigfiltz_normbaseline);

% Make signal processing trace figure
close all
injexampletraces = tiledlayout(2, 1, 'Padding','compact', 'TileSpacing','compact');
    
% Plot normalized signal trace - Saline 10
nexttile;
hold on;
plot(data(plotexampleSidx).sigfiltz_normbaseline, 'Color', S7color, 'LineWidth',0.3);
xline(data(plotexampleSidx).injt(1),'--','Color','#000000','FontSize',8,'LineWidth',1.5)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;
    
% Plot normalized signal trace - Morphine 10
nexttile;
hold on;
plot(data(plotexampleMidx).sigfiltz_normbaseline, 'Color', M8color, 'LineWidth',0.2);
xline(data(plotexampleMidx).injt(1),'--','Color','#000000','FontSize',8,'LineWidth',1.5)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;

% Set axis variables
allAxes = findall(gcf,'Type','axes');
set(allAxes, 'TickDir', 'out')
set(allAxes, 'XColor', '#000000', 'YColor', '#000000'); % Make axis lines black
set(gcf,'Color', 'w')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 5]);
exportgraphics(gcf,append(paperfigurepath,'Figure2_MorphineTraces_GCaMP6f','.svg'),'ContentType', 'vector','Resolution',300)


%% Figure 2 - VTA GCaMP6f Example Morphine Traces (Bin 24) with Transients
plotexampleSidx = 31; % 447 S 10
plotexampleMidx = 32; % 447 M 10

% Prep plot variables
currxstart = floor(data(plotexampleMidx).fs*60*69);
currxend = floor(data(plotexampleMidx).fs*60*72)-1;
currxlength = floor(data(plotexampleMidx).fs*180);    
currxticklabels = 69:1:72;
currxticks = floor((0:1:3).*60.*data(plotexampleMidx).fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz = max(data(plotexampleMidx).sigfiltz_normbaseline(currxstart:currxend));
yminsigfiltz = min(data(plotexampleMidx).sigfiltz_normbaseline(currxstart:currxend));

% Prep transients
plotexampleStransientlocs = transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc>currxstart) - currxstart;
plotexampleStransientvals = transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc>currxstart);

plotexampleMtransientlocs = transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc>currxstart) - currxstart;
plotexampleMtransientvals = transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc>currxstart);

% Make signal processing trace figure
close all
injtransexampletraces = tiledlayout(2, 1, 'Padding','compact', 'TileSpacing','compact');
    
% Plot normalized signal trace - Saline 10
nexttile;
hold on;
plot(data(plotexampleSidx).sigfiltz_normbaseline(currxstart:currxend), 'Color', S7color, 'LineWidth',0.4);
plot(plotexampleStransientlocs,plotexampleStransientvals,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;
    
% Plot normalized signal trace - Morphine 10
nexttile;
hold on;
plot(data(plotexampleMidx).sigfiltz_normbaseline(currxstart:currxend), 'Color', M8color, 'LineWidth',0.4);
plot(plotexampleMtransientlocs,plotexampleMtransientvals,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;

% Set axis variables
allAxes = findall(gcf,'Type','axes');
set(allAxes, 'TickDir', 'out')
set(allAxes, 'XColor', '#000000', 'YColor', '#000000'); % Make axis lines black
set(gcf,'Color', 'w')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 6.5, 5]);
exportgraphics(gcf,append(paperfigurepath,'Figure2_MorphineTracesTransientExample_GCaMP6f_Bin24','.svg'),'ContentType', 'vector','Resolution',300)


%% Figure 3 - NAcLS dLight1.3b Example Morphine Traces
plotexampleSidx = 155; % 447 S 10
plotexampleMidx = 156; % 434 M 10

% Prep plot variables
currxlength = length(data(plotexampleMidx).sigfiltz_normbaseline);    
currxticklabels = 0:6:72;
currxticks = floor(currxticklabels.*60.*data(plotexampleMidx).fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz = max(data(plotexampleMidx).sigfiltz_normbaseline);
yminsigfiltz = min(data(plotexampleMidx).sigfiltz_normbaseline);

% Make signal processing trace figure
close all
injexampletraces = tiledlayout(2, 1, 'Padding','compact', 'TileSpacing','compact');
    
% Plot normalized signal trace - Saline 10
nexttile;
hold on;
plot(data(plotexampleSidx).sigfiltz_normbaseline, 'Color', S7color, 'LineWidth',0.3);
xline(data(plotexampleSidx).injt(1),'--','Color','#000000','FontSize',8,'LineWidth',1.5)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;
    
% Plot normalized signal trace - Morphine 10
nexttile;
hold on;
plot(data(plotexampleMidx).sigfiltz_normbaseline, 'Color', M8color, 'LineWidth',0.2);
xline(data(plotexampleMidx).injt(1),'--','Color','#000000','FontSize',8,'LineWidth',1.5)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;

% Set axis variables
allAxes = findall(gcf,'Type','axes');
set(allAxes, 'TickDir', 'out')
set(allAxes, 'XColor', '#000000', 'YColor', '#000000'); % Make axis lines black
set(gcf,'Color', 'w')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 5]);
exportgraphics(gcf,append(paperfigurepath,'Figure3_MorphineTraces_dLight','.svg'),'ContentType', 'vector','Resolution',300)


%% Figure 3 - NAcLS dLight1.3b Example Morphine Traces (Bin 24) with Transients
plotexampleSidx = 155; % 447 S 10
plotexampleMidx = 156; % 434 M 10

% Prep plot variables
currxstart = floor(data(plotexampleMidx).fs*60*69);
currxend = floor(data(plotexampleMidx).fs*60*72)-1;
currxlength = floor(data(plotexampleMidx).fs*180);    
currxticklabels = 69:1:72;
currxticks = floor((0:1:3).*60.*data(plotexampleMidx).fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz = max(data(plotexampleMidx).sigfiltz_normbaseline(currxstart:currxend));
yminsigfiltz = min(data(plotexampleMidx).sigfiltz_normbaseline(currxstart:currxend));

% Prep transients
plotexampleStransientlocs = transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc>currxstart) - currxstart;
plotexampleStransientvals = transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc>currxstart);

plotexampleMtransientlocs = transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc>currxstart) - currxstart;
plotexampleMtransientvals = transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc>currxstart);

% Make signal processing trace figure
close all
injtransexampletraces = tiledlayout(2, 1, 'Padding','compact', 'TileSpacing','compact');
    
% Plot normalized signal trace - Saline 10
nexttile;
hold on;
plot(data(plotexampleSidx).sigfiltz_normbaseline(currxstart:currxend), 'Color', S7color, 'LineWidth',0.4);
plot(plotexampleStransientlocs,plotexampleStransientvals,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;
    
% Plot normalized signal trace - Morphine 10
nexttile;
hold on;
plot(data(plotexampleMidx).sigfiltz_normbaseline(currxstart:currxend), 'Color', M8color, 'LineWidth',0.4);
plot(plotexampleMtransientlocs,plotexampleMtransientvals,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;

% Set axis variables
allAxes = findall(gcf,'Type','axes');
set(allAxes, 'TickDir', 'out')
set(allAxes, 'XColor', '#000000', 'YColor', '#000000'); % Make axis lines black
set(gcf,'Color', 'w')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 6.5, 5]);
exportgraphics(gcf,append(paperfigurepath,'Figure3_MorphineTracesTransientExample_Bin24_dLight','.svg'),'ContentType', 'vector','Resolution',300)

%% Figure 4 - NAcLS GRABDA2h Example Morphine Traces
plotexampleSidx = 95; % 445 S 10
plotexampleMidx = 96; % 445 M 10

% Prep plot variables
currxlength = length(data(plotexampleMidx).sigfiltz_normbaseline);    
currxticklabels = 0:6:72;
currxticks = floor(currxticklabels.*60.*data(plotexampleMidx).fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz = max(data(plotexampleSidx).sigfiltz_normbaseline);
yminsigfiltz = min(data(plotexampleSidx).sigfiltz_normbaseline);

% Make signal processing trace figure
close all
injexampletraces = tiledlayout(2, 1, 'Padding','compact', 'TileSpacing','compact');
    
% Plot normalized signal trace - Saline 10
nexttile;
hold on;
plot(data(plotexampleSidx).sigfiltz_normbaseline, 'Color', S7color, 'LineWidth',0.3);
xline(data(plotexampleSidx).injt(1),'--','Color','#000000','FontSize',8,'LineWidth',1.5)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;
    
% Plot normalized signal trace - Morphine 10
nexttile;
hold on;
plot(data(plotexampleMidx).sigfiltz_normbaseline, 'Color', M8color, 'LineWidth',0.2);
xline(data(plotexampleMidx).injt(1),'--','Color','#000000','FontSize',8,'LineWidth',1.5)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;

% Set axis variables
allAxes = findall(gcf,'Type','axes');
set(allAxes, 'TickDir', 'out')
set(allAxes, 'XColor', '#000000', 'YColor', '#000000'); % Make axis lines black
set(gcf,'Color', 'w')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 5]);
exportgraphics(gcf,append(paperfigurepath,'Figure4_MorphineTraces_GRABDA2H','.svg'),'ContentType', 'vector','Resolution',300)


%% Figure 4 - NAcLS GRABDA2h Example Morphine Traces (Bin 24) with Transients
plotexampleSidx = 95; % 445 S 10
plotexampleMidx = 96; % 445 M 10

% Prep plot variables
currxstart = floor(data(plotexampleMidx).fs*60*69);
currxend = floor(data(plotexampleMidx).fs*60*72)-1;
currxlength = floor(data(plotexampleMidx).fs*180);    
currxticklabels = 69:1:72;
currxticks = floor((0:1:3).*60.*data(plotexampleMidx).fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz = max(data(plotexampleSidx).sigfiltz_normbaseline(currxstart:currxend));
yminsigfiltz = min(data(plotexampleSidx).sigfiltz_normbaseline(currxstart:currxend));

% Prep transients
plotexampleStransientlocs = transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc>currxstart) - currxstart;
plotexampleStransientvals = transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleSidx).transientquantification.maxloc>currxstart);

plotexampleMtransientlocs = transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc>currxstart) - currxstart;
plotexampleMtransientvals = transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleMidx).transientquantification.maxloc>currxstart);

% Make signal processing trace figure
close all
injtransexampletraces = tiledlayout(2, 1, 'Padding','compact', 'TileSpacing','compact');
    
% Plot normalized signal trace - Saline 10
nexttile;
hold on;
plot(data(plotexampleSidx).sigfiltz_normbaseline(currxstart:currxend), 'Color', S7color, 'LineWidth',0.4);
plot(plotexampleStransientlocs,plotexampleStransientvals,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;
    
% Plot normalized signal trace - Morphine 10
nexttile;
hold on;
plot(data(plotexampleMidx).sigfiltz_normbaseline(currxstart:currxend), 'Color', M8color, 'LineWidth',0.4);
plot(plotexampleMtransientlocs,plotexampleMtransientvals,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength]);
xticks(currxticks);
xticklabels(currxticklabels);
ylim([yminsigfiltz ymaxsigfiltz]);
xlabel('');
ylabel('');
hold off;

% Set axis variables
allAxes = findall(gcf,'Type','axes');
set(allAxes, 'TickDir', 'out')
set(allAxes, 'XColor', '#000000', 'YColor', '#000000'); % Make axis lines black
set(gcf,'Color', 'w')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 6.5, 5]);
exportgraphics(gcf,append(paperfigurepath,'Figure4_MorphineTracesTransientExample_Bin24_GRABDA2H','.svg'),'ContentType', 'vector','Resolution',300)

%% FIGURE 5 - Transient Color Grade Maps
plotdose = 10;

colormapyscalemins = [-5, -4, -5];
colormapyscalemaxes = [10, 6, 11];

colormapyscaleticks_GCaMP6f = [-5, 0, 5, 10];
colormapyscaleticks_GRABDA2h = [-3, 0, 3, 6];
colormapyscaleticks_dLight = [-5, 0, 5, 10];

colormapyscaleticks = [colormapyscaleticks_GCaMP6f; colormapyscaleticks_GRABDA2h; colormapyscaleticks_dLight];


for eachvirus = 1:length(viruses)
    currvirus = viruses{eachvirus};
    curryscalemin = colormapyscalemins(eachvirus);
    curryscalemax = colormapyscalemaxes(eachvirus);
    curryscaleticks = colormapyscaleticks(eachvirus,:);

    currfilepath = append(paperfigurepath,'Figure5_TransientColorMap_AllEvents_sigfiltznormBL_global_',currvirus,'_',num2str(plotdose));
    currfilepath_svg = append(paperfigurepath,'Figure5_TransientColorMap_AllEvents_sigfiltznormBL_global_',currvirus,'_',num2str(plotdose),'.',outputfiletype);
    curralltransientstruct = alltransientstruct_sigfiltnormzBL(strcmp({alltransientstruct_sigfiltnormzBL.Virus},currvirus) & strcmp([alltransientstruct_sigfiltnormzBL.Phase],'POSTINJ') ...
        & [alltransientstruct_sigfiltnormzBL.Dose] == plotdose);
    
    close all   
    alltransientcolormaps = plotTransientColorGradeGroups(curralltransientstruct,'InjType','transientstreamdata',customColorScale,currvirus,'climtype','global', 'yscaletype', 'manual', 'yscalemin', curryscalemin, 'yscalemax', curryscalemax);

    nexttile(1)
    colorbar off
    title('')
    yticks(curryscaleticks)
    
    nexttile(2)
    set(gca, 'YTickLabel', [])    
    ylabel('')
    title('')
    yticks(curryscaleticks)

    cb = colorbar;
    cb.Box = 'off';
    
    title(alltransientcolormaps, '')
    allaxes = findall(gcf,'Type','axes');
    set(allaxes, 'XColor', '#000000', 'YColor', '#000000', 'Color', 'none', 'TickDir', 'out', 'LineWidth', .85, 'FontSize', 7, 'TickLength', [0.03 0.03])
    
    set(gcf, 'Units', 'inches', 'Position', [0, 0, 2.9, 1.8]);
    exportgraphics(gcf,append(currfilepath,'.','eps'),'Resolution',300)

end

%% Figure 1 - Example Saline Traces with Transients
plotexampleSidx_GCaMP6f = find([data.SubjectID]==429 & strcmp({data.InjType},'S')&[data.Dose]==10); % 429 S 10
plotexampleSidx_dLight =  find([data.SubjectID]==479 & strcmp({data.InjType},'S')&[data.Dose]==10); % 479 S 10
plotexampleSidx_GRABDA2h = find([data.SubjectID]==449 & strcmp({data.InjType},'S')&[data.Dose]==10); % 449 S 5

fs = data(1).fs;
tracesecs = 30;


% Prep plot variables - GCaMP6f
currxstart_GCaMP6f = floor(fs*4221);
currxend_GCaMP6f = currxstart_GCaMP6f+floor(tracesecs*fs);
currxlength_GCaMP6f = floor(fs*tracesecs);    
currxticklabels_GCaMP6f = 0:5:tracesecs;
currxticks_GCaMP6f = floor(currxticklabels_GCaMP6f*fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz_GCaMP6f = max(data(plotexampleSidx_GCaMP6f).sigfiltz_normbaseline(currxstart_GCaMP6f:currxend_GCaMP6f));
yminsigfiltz_GCaMP6f = min(data(plotexampleSidx_GCaMP6f).sigfiltz_normbaseline(currxstart_GCaMP6f:currxend_GCaMP6f));

plotexampleStransientlocs_GCaMP6f = transientdata_sigfiltznormBL(plotexampleSidx_GCaMP6f).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleSidx_GCaMP6f).transientquantification.maxloc>currxstart_GCaMP6f) - currxstart_GCaMP6f;
plotexampleStransientvals_GCaMP6f = transientdata_sigfiltznormBL(plotexampleSidx_GCaMP6f).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleSidx_GCaMP6f).transientquantification.maxloc>currxstart_GCaMP6f);

% Prep plot variables - dLight
currxstart_dLight = floor(fs*3553);
currxend_dLight = currxstart_dLight+floor(tracesecs*fs);
currxlength_dLight = floor(fs*tracesecs);    
currxticklabels_dLight = 0:5:tracesecs;
currxticks_dLight = floor(currxticklabels_dLight*fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz_dLight = max(data(plotexampleSidx_dLight).sigfiltz_normbaseline(currxstart_dLight:currxend_dLight));
yminsigfiltz_dLight = min(data(plotexampleSidx_dLight).sigfiltz_normbaseline(currxstart_dLight:currxend_dLight));

plotexampleStransientlocs_dLight = transientdata_sigfiltznormBL(plotexampleSidx_dLight).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleSidx_dLight).transientquantification.maxloc>currxstart_dLight) - currxstart_dLight;
plotexampleStransientvals_dLight = transientdata_sigfiltznormBL(plotexampleSidx_dLight).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleSidx_dLight).transientquantification.maxloc>currxstart_dLight);

% Prep plot variables - GRABDA2h
currxstart_GRABDA2h = floor(fs*717);
currxend_GRABDA2h = currxstart_GRABDA2h+floor(tracesecs*fs);
currxlength_GRABDA2h = floor(fs*tracesecs);    
currxticklabels_GRABDA2h = 0:5:tracesecs;
currxticks_GRABDA2h = floor(currxticklabels_GRABDA2h*fs); % Determine x axis ticks - add ticks every 5 minutes

ymaxsigfiltz_GRABDA2h = max(data(plotexampleSidx_GRABDA2h).sigfiltz_normbaseline(currxstart_GRABDA2h:currxend_GRABDA2h));
yminsigfiltz_GRABDA2h = min(data(plotexampleSidx_GRABDA2h).sigfiltz_normbaseline(currxstart_GRABDA2h:currxend_GRABDA2h));

plotexampleStransientlocs_GRABDA2h = transientdata_sigfiltznormBL(plotexampleSidx_GRABDA2h).transientquantification.maxloc(transientdata_sigfiltznormBL(plotexampleSidx_GRABDA2h).transientquantification.maxloc>currxstart_GRABDA2h) - currxstart_GRABDA2h;
plotexampleStransientvals_GRABDA2h = transientdata_sigfiltznormBL(plotexampleSidx_GRABDA2h).transientquantification.maxval(transientdata_sigfiltznormBL(plotexampleSidx_GRABDA2h).transientquantification.maxloc>currxstart_GRABDA2h);

% Make signal processing trace figure
close all
injtransexampletraces = tiledlayout(3, 1, 'Padding','compact', 'TileSpacing','compact');
    
% Plot GCaMP6f
nexttile;
hold on;
plot(data(plotexampleSidx_GCaMP6f).sigfiltz_normbaseline(currxstart_GCaMP6f:currxend_GCaMP6f), 'Color', '#008121', 'LineWidth',0.4);
plot(plotexampleStransientlocs_GCaMP6f,plotexampleStransientvals_GCaMP6f,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength_GCaMP6f]);
xticks(currxticks_GCaMP6f);
xticklabels(currxticklabels_GCaMP6f);
ylim([yminsigfiltz_GCaMP6f ymaxsigfiltz_GCaMP6f]);
xlabel('');
ylabel('');
hold off;
    
% Plot dLight
nexttile;
hold on;
plot(data(plotexampleSidx_dLight).sigfiltz_normbaseline(currxstart_dLight:currxend_dLight), 'Color', '#003FD1', 'LineWidth',0.4);
plot(plotexampleStransientlocs_dLight,plotexampleStransientvals_dLight,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength_dLight]);
xticks(currxticks_dLight);
xticklabels(currxticklabels_dLight);
ylim([yminsigfiltz_dLight ymaxsigfiltz_dLight]);
xlabel('');
ylabel('');
hold off;

% Plot GRABDA2h
nexttile;
hold on;
plot(data(plotexampleSidx_GRABDA2h).sigfiltz_normbaseline(currxstart_GRABDA2h:currxend_GRABDA2h), 'Color', '#EA007D', 'LineWidth',0.4);
plot(plotexampleStransientlocs_GRABDA2h,plotexampleStransientvals_GRABDA2h,'o', 'Color', '#000000','MarkerSize', 3, 'LineWidth', .6)
xlim([0 currxlength_GRABDA2h]);
xticks(currxticks_GRABDA2h);
xticklabels(currxticklabels_GRABDA2h);
ylim([yminsigfiltz_GRABDA2h ymaxsigfiltz_GRABDA2h]);
xlabel('');
ylabel('');
hold off;

% Set axis variables
allAxes = findall(gcf,'Type','axes');
set(allAxes, 'TickDir', 'out')
set(allAxes, 'XColor', '#000000', 'YColor', '#000000'); % Make axis lines black
set(gcf,'Color', 'w')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 5, 7.3]);
exportgraphics(gcf,append(paperfigurepath,'SupplementaryFigure_Control_SensorTransientExampleTraces','.svg'),'ContentType', 'vector','Resolution',300)

