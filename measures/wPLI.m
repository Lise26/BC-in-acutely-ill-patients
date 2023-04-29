% wPLI class (herits from associationMeasure)

% Code coming from Ghita Ait Ouhmane's Master Thesis
% Ait Ouhmane, Ghita. 2022. “Assessing consciousness at the Intensive Care
% Unit (ICU)." Master s thesis, Université Libre de Bruxelles.

classdef wPLI < frequencyDomain
    methods
       % Constructor
       function obj = wPLI(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize)
            obj@frequencyDomain(windowLength, windowOverlap,fhband,flband, nbIndicators, outSize);
        end

        function wpli = measure(obj, x, y, fs)
            % wPLI as defined by Vinck et al.
            % For each epoch, a Hann window function is used on the signal 
            % and the wPLI is computed using cpsd 

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

            % The cpsd is averaged across epochs which gives a wPLI array 
            % for frequencies in [0,fs/2]
            num = abs(mean(imag(windowValues),1));
            den = mean(abs(imag(windowValues)),1);
            wpliVinck = num./den;

            wpli = zeros(obj.freqs,1);
            for b = 1:obj.freqs
                [range,~] = find(fw(:,1)>=obj.fhbands(b) & fw(:,1)<obj.flbands(b));
                wpli(b) = mean(wpliVinck(range),'omitnan');
            end
       end

   end
end
