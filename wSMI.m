% weighted Symboli Mutual Information class (herits from timeDomain)

classdef wSMI < timeDomain
    properties
        tau {mustBeNumeric}
        kernel {mustBeNumeric}
        symbols {mustBeNumeric}
    end

    methods
       % Constructor
       function obj = wSMI(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize, tau, kernel)
            obj@timeDomain(windowLength, windowOverlap,fhband,flband,nbIndicators,outSize);
            obj.tau = tau;
            obj.kernel = kernel;
            obj.symbols = factorial(kernel);
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

        function wsmi = association(obj,x,y)
            % Symbolic transform
            xHat = obj.symbolicTransform(x);
            yHat = obj.symbolicTransform(y);
            
            % weighted Mutual information
            hxy = histcounts2(xHat, yHat, obj.symbols);
            hx = histcounts(xHat, obj.symbols);
            hy = histcounts(yHat, obj.symbols);
            pxy = hxy/sum(hxy,"all");
            px = hx/sum(hx);
            py = hy/sum(hy);
            
            weights = [0,1,1,1,0,1;1,0,1,1,1,0;1,1,0,0,1,1;1,1,0,0,1,1;0,1,1,1,0,1;1,0,1,1,1,0];
            wsmi = 0;
            for i = 1:obj.symbols
                for j = 1:obj.symbols
                    if pxy(i,j) > 0
                        wsmi = wsmi + weights(i,j)*pxy(i,j)*log10(pxy(i,j)/(px(i)*py(j)));
                    end
                end
            end
            wsmi = (1/log10(obj.symbols)) * wsmi;
        end

        function symbolic = symbolicTransform(obj,signal)
            symbolic = zeros(length(signal)-obj.tau*(obj.kernel-1),1);
            for s = 1:length(symbolic)
                range = s:obj.tau:s+obj.kernel*obj.tau;
                if (signal(range(2))<signal(range(1))) && (signal(range(2))>signal(range(3))) % negative line
                    symbolic(s) = 1;
                elseif (signal(range(3))<signal(range(1))) && (signal(range(3))>signal(range(2))) % U left
                    symbolic(s) = 2;
                elseif (signal(range(1))<signal(range(2))) && (signal(range(1))>signal(range(3))) % bridge right
                    symbolic(s) = 3;
                elseif (signal(range(1))<signal(range(3))) && (signal(range(1))>signal(range(2))) % U right
                    symbolic(s) = 4;
                elseif (signal(range(2))<signal(range(3))) && (signal(range(2))>signal(range(1))) % positive line
                    symbolic(s) = 5;
                elseif (signal(range(3))<signal(range(2))) && (signal(range(3))>signal(range(1))) % bridge left
                    symbolic(s) = 6;
                end
            end
        end

        function res = processEpochs(obj, windowValues)
            res = mean(windowValues);
        end

   end
end
