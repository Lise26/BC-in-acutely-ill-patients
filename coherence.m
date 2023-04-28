% coherence class (herits from associationMeasure)
% Computes the coherence

classdef coherence < frequencyDomain
    methods
        % Constructor
        function obj = coherence(windowLength, windowOverlap, fhband, flband, nbIndicators)
            outSize = 129;
            obj@frequencyDomain(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize);
        end

        function res = measure(obj,x,y,fs)
            nbWindows = floor(length(x)/obj.windowStep)-2;
            windowValues = zeros(nbWindows,obj.outSize);
            fValues = zeros(obj.outSize,nbWindows);
            parfor w = 1:nbWindows
                start = obj.windowStep*w;
                % Windowing
                xw = x(start:start+obj.windowLength,1);
                yw = y(start:start+obj.windowLength,1);
                % Normalization
                xw = (xw - mean(xw))/std(xw);
                yw = (yw - mean(yw))/std(yw);
                [windowValues(w,:), fValues(:,w)] = obj.association(xw,yw,fs);
            end
            res = obj.processEpochs(windowValues,fValues);
        end

        % Coherency computation
        function [coherency,fxy] = association(obj,x,y,fs)
            xw = hann(obj.windowLength).'.*x(1:end-1).';
            yw = hann(obj.windowLength).'.*y(1:end-1).';
            [sxy, fxy] = cpsd(xw,yw,[],[],[],fs);
            [sxx, ~] = cpsd(xw,xw,[],[],[],fs);
            [syy, ~] = cpsd(yw,yw,[],[],[],fs);

            coherency = sxy./sqrt(sxx.*syy);
        end

        function res = processEpochs(obj, windowValues, f)
            res = zeros(obj.nbIndicators,obj.freqs);
            meanWind = mean(windowValues,1);

            coherence = abs(meanWind);
            imagCoh = imag(meanWind);
            for b = 1:obj.freqs
                [range,~] = find(f(:,1)>=obj.fhbands(b) & f(:,1)<obj.flbands(b));
                % Mean coherence
                res(1,b) = mean(coherence(range),'omitnan');
                % Max coherence
                res(2,b) = max(abs(coherence(range)));
                % Imaginary part of the coherency
                res(3,b) = mean(imagCoh(range), 'omitnan');
            end

        end
    end
end
