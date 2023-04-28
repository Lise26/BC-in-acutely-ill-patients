% EEGDATA CLASS

% Manages the raw EEG data
%
% Functions to extract the relevant information from the EEG data Matlab
% structure and to apply pre-processing steps for the removal of the 
% unwanted components possibly present in the signals (noise, artifacts, 
% interferences and background activity).

classdef EEGData
    properties
        Channels {mustBeNumeric}
        Points {mustBeNumeric}
        SampleFrequency {mustBeNumeric}
        Reference {mustBeNumeric}
        ChannelsIdx
        Data
        FreqData
    end

    methods
        % Constructor
        function obj = EEGData(EEGFile, startTimePoint, timeToAnalyze)
            load(EEGFile, "EEG");
            if strcmp(timeToAnalyze, "all")
                obj.Points = size(EEG.data,1);
            else
                obj.Points = timeToAnalyze;
            end
            EEGChannels = ["F3","C3","P3","Cz","F4","C4","P4","Fp1","Fp2","F7","T3","T5","O1","O2","T6","T4","F8","Fz","Pz","A1","A2"];
            for i = 1:size(EEG.data, 2)
                if any(strcmp(string(EEG.channels(i).Name), EEGChannels))
                    obj.ChannelsIdx = [obj.ChannelsIdx, string(EEG.channels(i).Name)];
                end
            end
            obj.Channels = length(obj.ChannelsIdx) - 2;
            obj.Data = EEG.data(startTimePoint:startTimePoint+obj.Points-1,1:obj.Channels+2);
            obj.SampleFrequency = EEG.channels.Fs;
        end

        % Re-referencing
        function obj = rereferenceToAverage(obj)
            data = obj.Data;
            data(:,obj.ChannelsIdx=="Fp1") = [];
            data(:,obj.ChannelsIdx=="Fp2") = [];
            data(:,obj.ChannelsIdx=="Cz") = [];
            average = mean(data,2);
            for i=1:obj.Channels
                obj.Data(:,i) = obj.Data(:,i) - average;
            end
        end

        function obj = rereferenceToMastoids(obj)
            average = mean(obj.Data(:,20:21),2);
            for i=1:obj.Channels
                obj.Data(:,i) = obj.Data(:,i) - average;
            end
            obj.Data(:,20:21) = [];
        end

        function obj = rereferenceToLaplacian(obj)
            laplace = {["F3","Fp1","F7","C3","Fz"],["C3","F3","T3","P3","Cz"],["P3","C3","T5","O1","Pz"],["F4","Fp2","Fz","C4","F8"],...
                ["Cz","Fz","C3","Pz","C4","F3","F4","P4","P3"],["C4","F4","Cz","P4","T4"],["P4","C4","Pz","O2","T6"],...
                ["Fp1","F7","F3","Fz","Fp2"],["Fp2","Fz","F4","F8","Fp1"],["F7","Fp1","F3","C3","T3"], ["T3","F7","C3","T5","A1"],...
                ["T5","T3","P3","O1","C3"],["O1","T5","P3","Pz","O2"],["O2","T6","P4","Pz","O1"],["T6","T4","P4","O2","C4"],...
                ["T4","F8","C4","T6","A2"],["F8","Fp2","F4","T4","C4"],["Fz","Fp1","Fp2","F3","Cz","F4"],["Pz","O1","O2","P3","Cz","P4"]};
            for i = 1:length(laplace)
                surroundings = length(laplace{i});
                if laplace{i}(1) == "Fz" || laplace{i}(1) == "Pz"
                    sum = (obj.Data(:,obj.ChannelsIdx==laplace{i}(2)) + obj.Data(:,obj.ChannelsIdx==laplace{i}(3)))/2;
                    for j = 4:surroundings
                        sum = sum + obj.Data(:,obj.ChannelsIdx==laplace{i}(j));
                    end
                    newReference = (1/(surroundings-2)) * sum;
                else
                    sum = 0;
                    for j = 2:surroundings
                        sum = sum + obj.Data(:,obj.ChannelsIdx==laplace{i}(j));
                    end
                    newReference = (1/(surroundings-1)) * sum;
                end
                obj.Data(:,obj.ChannelsIdx==laplace{i}(1)) = obj.Data(:,obj.ChannelsIdx==laplace{i}(1)) - newReference;
            end
        end

        % Filtering
        function obj = filtering(obj, lowFreqs, highFreqs)
            for band = 1:length(lowFreqs)
                [bHigh, aHigh] = butter(1,lowFreqs(band)/(obj.SampleFrequency/2),"high");    % 1st order highpass
                [bLow, aLow] = butter(2,highFreqs(band)/(obj.SampleFrequency/2));            % 2nd order lowpass
                for i=1:obj.Channels
                    filtHData = filtfilt(bHigh,aHigh,obj.Data(:,i));
                    obj.FreqData(:,i,band) = filtfilt(bLow,aLow,filtHData);
                end
            end
        end

    end
end