% correctedCrossCorrelation class (herits from timeDomain)
% Computes the corrected cross-correlation

classdef correctedCrossCorrelation < timeDomain
    properties
        maxLag {mustBeNumeric}
    end

    methods
        % Constructor
        function obj = correctedCrossCorrelation(windowLength, windowOverlap, fhband, flband, nbIndicators, maxLag)
            outSize = maxLag;
            obj@timeDomain(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize);
            obj.maxLag = maxLag;
        end

        function res = measure(obj, x, y)
            nbWindows = floor(length(x)/obj.windowStep)-2;
            windowValues = zeros(nbWindows-2,obj.outSize);
            parfor w = 2:nbWindows-1
                start = obj.windowStep*w;
                xw = x(start:start+obj.windowLength,1);
                yw = y(start-obj.maxLag:start+obj.windowLength+obj.maxLag,1);
                windowValues(w-1,:) = obj.association(xw,yw)';
            end
            res = obj.processEpochs(windowValues);
        end

        % Corrected cross-correlation
        function corr = association(obj,x,y)
            cross = crossCorrelation(obj.windowLength,obj.windowOverlap,obj.fhbands,obj.flbands,obj.nbIndicators,obj.maxLag);
            cc = cross.association(x,y);
            corr = zeros(obj.outSize,1);
            for lag=1:obj.maxLag
                corr(lag,1) = (1/2) * (cc(lag+obj.maxLag+1) - cc(-lag+obj.maxLag+1));
            end
        end

        function res = processEpochs(obj,windowValues)
            meanWind = mean(windowValues,1);
            [res(1), res(3)] = max(abs(meanWind));
            res(2) = mean(meanWind);
        end

    end
end