% associationMeasure class

% Parent class of all association measures

classdef associationMeasure
    properties
        windowLength {mustBeNumeric}
        windowOverlap {mustBeNumeric}
        windowStep {mustBeNumeric}
        fhbands
        flbands
        freqs {mustBeNumeric}
        nbIndicators {mustBeNumeric}
        outSize {mustBeNumeric}
    end

    methods
        % Constructor
        function obj = associationMeasure(windowLength, windowOverlap, fhband, flband, nbIndicators, outSize)
            obj.windowLength = windowLength;
            obj.windowOverlap = windowOverlap;
            obj.windowStep = windowLength - windowOverlap;
            obj.fhbands = fhband;
            obj.flbands = flband;
            obj.freqs = length(fhband);
            obj.nbIndicators = nbIndicators;
            obj.outSize = outSize;
        end

    end
end
