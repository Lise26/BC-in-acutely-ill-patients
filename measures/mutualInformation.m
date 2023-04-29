% mutualInformation class (herits from timeDomain)

classdef mutualInformation < timeDomain
    properties
        bins {mustBeNumeric}
    end

    methods
       % Constructor
       function obj = mutualInformation(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize, bins)
            obj@timeDomain(windowLength, windowOverlap,fhband,flband,nbIndicators, outSize);
            obj.bins = bins;
       end

       function res = measure(obj, x, y)
            nbWindows = floor(length(x)/obj.windowStep)-2;
            windowValues = zeros(nbWindows,obj.outSize);
            parfor w = 1:nbWindows
                start = obj.windowStep*w;
                % Windowing
                xw = x(start:start+obj.windowLength,1);
                yw = y(start:start+obj.windowLength,1);
                % Normalization
                xw = (xw - mean(xw))/std(xw);
                yw = (yw - mean(yw))/std(yw);
                windowValues(w) = obj.association(xw,yw);
            end
            res = obj.processEpochs(windowValues);
        end

        function mi = association(obj,x,y)
            hxy = histcounts2(x, y, obj.bins);
            hx = histcounts(x, obj.bins);
            hy = histcounts(y, obj.bins);
            
            pxy = hxy/sum(hxy,"all");
            px = hx/sum(hx);
            py = hy/sum(hy);
            
            mi = 0;
            for i = 1:obj.bins
                for j = 1:obj.bins
                    if pxy(i,j) > 0
                        mi = mi + pxy(i,j)*log10(pxy(i,j)/(px(i)*py(j)));
                    end
                end
            end
        end

        function res = processEpochs(obj, windowValues)
            res = mean(windowValues);
        end

   end
end
