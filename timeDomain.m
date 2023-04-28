% timeDomain class (herits from associationMeasure)

classdef timeDomain < associationMeasure
    methods
        function obj = timeDomain(windowLength, windowOverlap,fhband, flband, nbIndicators, outSize)
            obj@associationMeasure(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize)
        end

        function mat = matrix(obj, data)
            mat = zeros(data.Channels, data.Channels, obj.nbIndicators, obj.freqs);
            data = data.filtering(obj.fhbands,obj.flbands);
            for b = 1:obj.freqs
                for i = 1:data.Channels
                    for j = i+1:data.Channels
                        mat(i,j,:,b) = obj.measure(data.FreqData(:,i,b), data.FreqData(:,j,b));
                    end
                end
            end
        end

    end
end
