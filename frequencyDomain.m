% frequencyDomain class (herits from associationMeasure)

classdef frequencyDomain < associationMeasure
    methods
        function obj = frequencyDomain(windowLength, windowOverlap,fhband,flband, nbIndicators, outSize)
            obj@associationMeasure(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize)
        end

        function mat = matrix(obj, data)
            mat = zeros(data.Channels, data.Channels, obj.nbIndicators, obj.freqs);
            data = data.filtering(obj.fhbands(end),obj.flbands(end));
            for i = 1:data.Channels
                for j = i+1:data.Channels
                    mat(i,j,:,:) = obj.measure(data.FreqData(:,i), data.FreqData(:,j), data.SampleFrequency);
                end
            end
        end
    end
end