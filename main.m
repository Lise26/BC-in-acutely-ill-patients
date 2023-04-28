% MAIN
clear; clc;

% GENERAL PARAMETERS
% The data of each patient is stored in a MAT structure with 2 fields:
% - channels: 1x21 structures with 2 fields: Name and Fs
% - data: 225250x21 double with the EEG voltages (in microVolts)
dataFolder = pwd + "\USI\DataUSI\Cohort1\MAT\";
nbPatients = 100;
startTime = 1;
timeToAnalyze = 300000;         % 20 minutes
fhband = [1,4,8,13,1,1];        % low-frequencies
flband = [4,8,13,20,13,20];     % high-frequencies 
band = ["Delta", "Theta", "Alpha", "Beta", "Quasi-Broad", "Broad"];
windowSize = 250;               % 1 second
overlap = 125;                  % 0.5 second
maxLag = 50;
bins = 256;

% MEASURES
indicators = 3;
% outName = "connectivityCC.mat";
% measure = crossCorrelation(windowSize,overlap,fhband,flband,indicators,maxLag);
% connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)
% 
outName = "connectivityCORR.mat";
maxlag = 100;
measure = correctedCrossCorrelation(windowSize,overlap,fhband,flband,indicators,maxLag);
connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)
% 
% outName = "connectivityCOH.mat";
% measure = coherence(windowSize,overlap,fhband,flband,indicators);
% connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)
% 
% indicators = 1;
% outSize = 129;
% outName = "connectivityPLI.mat";
% measure = PLI(windowSize, overlap,fhband,flband,indicators, outSize);
% connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)

% indicators = 1;
% outSize = 129;
% outName = "connectivityPLV.mat";
% measure = PLV(windowSize, overlap,fhband,flband,indicators, outSize);
% connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)

% outSize = 129;
% outName = "connectivityWPLI.mat";
% measure = wPLI(windowSize, overlap,fhband,flband,indicators, outSize);
% connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)

% outSize = 1;
% bins = 256;
% outName = "connectivityMI.mat";
% measure = mutualInformation(windowSize,overlap,fhband,flband,indicators,outSize,bins);
% connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)
% 
% tau = 4;
% kernel = 3;
% outName = "connectivityWSMI.mat";
% measure = wSMI(windowSize,overlap,fhband,flband,indicators,outSize,tau,kernel);
% connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pipeline to build connectivity networks
%
% Starts from EEG data to end up with a connectivity indicator of the brain
% network by averaging the connectivity matrix
% 
% STEPS:
%   1: Data preprocessing
%       - re-referencing
%       - filtering
%   2: Network computation [in every frequency band]
%       - windows processing
%       - normalization
%       - associations computation
%   3: Average matrix per frequency band

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function connectivityMatrix(dataFolder,nbPatients,startTime,timeToAnalyze,band,outName,measure)
    connectivity = struct();

    for p = 1:nbPatients
        disp(" ----- Processing patient n " + p)
        tstart = datetime('now','Format','HH:mm:ss.SSS');
        
        fileName = dataFolder + p + ".mat";
        data = EEGData(fileName, startTime, timeToAnalyze);
    
        % 1: PREPROCESSING
        % Re-referencing to mastoids
        data = data.rereferenceToMastoids;
    
        % 2: NETWORK COMPUTATION
        assocMatrix = measure.matrix(data);
        for ind = 1:measure.nbIndicators
            for b = 1:measure.freqs
                brainNetwork = network(assocMatrix(:,:,ind,b));
                connectivity(ind).patient(p).matrix(:,:,b) = brainNetwork.Edges;
        
                % Displays the matrices per frequency bands
                % figure
                % heatmap(data.ChannelsIdx(1:19),data.ChannelsIdx(1:19),brainNetwork.Edges,'ColorMap',jet,'ColorLimits',[0 1]);
                % title("wPLI - Band: "+ band(b))
                % set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
                % 3: AVERAGE
                for i = 1:19
                    brainNetwork.Edges(i,i,b) = nan;
                end
                connectivity(ind).patient(p).average(b) = mean(brainNetwork.Edges, 'all', 'omitnan');
            end
        end
        tend = datetime('now','Format','HH:mm:ss.SSS');
        ttot = tend- tstart; % total processing time
        ttot
    
    end
    save(outName, "connectivity");
end

