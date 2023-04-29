% crossCorrelation class (herits from timeDomain)
% Computes the cross-correlation

classdef crossCorrelation < timeDomain
    properties
        maxLag {mustBeNumeric}
    end

    methods
        % Constructor
        function obj = crossCorrelation(windowLength, windowOverlap, fhband, flband, nbIndicators, maxLag)
            outSize = 2*maxLag + 1;
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

        % Cross-correlation computation
        function cc = association(obj,x,y)
            cc = zeros(obj.outSize,1);
            meanx = mean(x);
            for lag = 1:2*obj.maxLag+1
                yc = y(lag:lag+obj.windowLength);
                meany = mean(yc);
                frac = 1/(std(x)*std(yc)*obj.windowLength);
                cross = 0;
                for t=1:obj.windowLength
                    cross = cross + ((x(t)-meanx)*(yc(t)-meany));
                end
                cc(lag,1) = frac*cross;
            end
        end

        function res = processEpochs(obj,windowValues)
            meanWind = mean(windowValues,1);
            [res(1), index] = max(abs(meanWind));
            res(2) = mean(meanWind);
            res(3) = index - obj.maxLag - 1;
        end

    end
end
