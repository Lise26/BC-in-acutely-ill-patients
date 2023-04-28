% PLV class (herits from associationMeasure)

classdef PLV < timeDomain
    methods
       % Constructor
       function obj = PLV(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize)
            obj@timeDomain(windowLength, windowOverlap,fhband,flband, nbIndicators, outSize);
        end

        function plv = measure(obj, x, y)
            nbWindows = floor(length(x)/obj.windowStep)-2;
            windowValues = zeros(nbWindows,1);
            parfor w = 1:nbWindows
                start = obj.windowStep*(w-1)+1;
                % Normalization
                xw = (x(start:start+obj.windowLength,1) - mean(x(start:start+obj.windowLength,1)))/std(x(start:start+obj.windowLength,1));
                yw = (y(start:start+obj.windowLength,1) - mean(y(start:start+obj.windowLength,1)))/std(y(start:start+obj.windowLength,1));
                xw = hann(obj.windowLength).'.*xw(1:end-1).';
                yw = hann(obj.windowLength).'.*yw(1:end-1).';
                
                xh = hilbert(xw);
                yh = hilbert(yw);
                for t = 1:obj.windowLength
                    phaseX = atan(imag(xh(t))/real(xh(t)));
                    phaseY = atan(imag(yh(t))/real(yh(t)));
                    phaseDiff = phaseY - phaseX;
                    windowValues(w) = windowValues(w) + exp(1i*phaseDiff);
                end
                windowValues(w) = (1/obj.windowLength)*abs(windowValues(w));
            end
            plv = mean(windowValues);
       end

   end
end
