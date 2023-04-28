% PLI class (herits from associationMeasure)

classdef PLI < frequencyDomain
    methods
       % Constructor
       function obj = PLI(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize)
            obj@frequencyDomain(windowLength, windowOverlap,fhband,flband, nbIndicators, outSize);
        end

        function pli = measure(obj, x, y, fs)
            nbWindows = floor(length(x)/obj.windowStep)-2;
            windowValues = zeros(nbWindows,obj.outSize);
            fw = zeros(obj.outSize, nbWindows);
            parfor w = 1:nbWindows
                start = obj.windowStep*(w-1)+1;
                % Normalization
                xw = (x(start:start+obj.windowLength,1) - mean(x(start:start+obj.windowLength,1)))/std(x(start:start+obj.windowLength,1));
                yw = (y(start:start+obj.windowLength,1) - mean(y(start:start+obj.windowLength,1)))/std(y(start:start+obj.windowLength,1));
    
                xw = hann(obj.windowLength).'.*xw(1:end-1).';
                yw = hann(obj.windowLength).'.*yw(1:end-1).';
                [windowValues(w,:), fw(:,w)] = cpsd(xw,yw,[],[],[],fs);
            end

            pliFreq = abs(mean(sign(imag(windowValues)),1));
            pli = zeros(obj.freqs,1);
            for b = 1:obj.freqs
                [range,~] = find(fw(:,1)>=obj.fhbands(b) & fw(:,1)<obj.flbands(b));
                pli(b) = mean(pliFreq(range),'omitnan');
            end
       end

   end
end
