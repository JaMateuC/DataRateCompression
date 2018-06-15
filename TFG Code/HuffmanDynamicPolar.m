function [error,avglen,signalSize] = HuffmanDynamicPolar(signal,numValuesAng,numValuesRad,trueValue,plots,huffman)

maxR = .65;
maxA = 6;
startSig = zeros(length(signal),1);

tmwaveform2 = normalization(signal);
startSig(1) = tmwaveform2(1);

distancePolar = abs(abs(tmwaveform2(1:end-1))-abs(tmwaveform2(2:end)));
distanceAngle = atan2(imag(tmwaveform2(1:end-1)),real(tmwaveform2(1:end-1))) - atan2(imag(tmwaveform2(2:end)),real(tmwaveform2(2:end)));

TreuValuesVecIndx = 1:trueValue:length(tmwaveform2)-trueValue;


intervalRad = maxR/(numValuesRad);
intervalRadVector = [0+intervalRad:intervalRad:maxR]';
intervalAng = 2*maxA/(numValuesAng-1);
intervalAngVector = [-maxA:intervalAng:maxA]';

minDAll = abs(distancePolar-intervalRadVector');
[~,minInd] = min(minDAll,[],2);
tmwaveformCRad = intervalRadVector(minInd)'.';

minDAll = abs(distanceAngle-intervalAngVector');
[~,minInd] = min(minDAll,[],2);
tmwaveformCAng = intervalAngVector(minInd)'.';

for i=1:length(distanceAngle)
    
    if(sum(TreuValuesVecIndx == (i+1)))
        startSig(i+1) = tmwaveform2(i+1);
    else
        SAng = atan2(imag(startSig(i)),real(startSig(i))) - tmwaveformCAng(i);
        SRad = abs(startSig(i)) + tmwaveformCRad(i);
        startSig(i+1) = SRad * cos(SAng) + 1i * SRad * sin(SAng);
    end
    
end

tmwaveformC = round(startSig,5);

error = EVM(tmwaveform2,tmwaveformC,plots);

if(plots)
    
    figure
    histogram(distancePolar,100)
    
    figure
    histogram(distanceAngle,100)
    
    intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
    figure
    plot(intervalVector(:,1),intervalVector(:,2), 'xr')
    xlabel('Phase')
    ylabel('Quadrature')
    title('Constellation')
    figure
    voronoi(intervalVector(:,1),intervalVector(:,2))
    title('Constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    
    figure
    plot(intervalAmplitudeVector,'x')
    title('Levels Assigments')
    ylabel('Interval')
    
    figure
    plot(tmwaveformC, 'x')
    title('Signal compressed')
    ylabel('Quadrature')
    xlabel('Phase')
    
    intervalVectorTog = intervalVector(:,1) + 1i * intervalVector(:,2);
    
    vectorSectionSize = ones(numValues,1)*(numValues);    
    
    minDAll = abs(tmwaveform2-intervalVectorTog');
    [~,minInd] = min(minDAll,[],2);
    
    errorsVector = [abs(tmwaveformC - tmwaveform2),minInd];
    errorsHist = zeros(length(intervalVectorTog),1);
    for i=1:length(intervalVectorTog)
        errorsHist(i) = sum(errorsVector(errorsVector(:,2) == i));
    end
    errorsHistCell = mat2cell(errorsHist',1,vectorSectionSize');
    errorsHistTotal = cellfun(@(x) sum(x)/sum(errorsHist)*100,errorsHistCell);
    figure
    bar(errorsHistTotal)
    title('Error acumulation on each level')
    ylabel('Error acumulated')
    xlabel('Interval')
    
end

if(huffman)
    intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
    intervalVector = round(intervalVector(1:numValues,2),5);
    AccSamp = histcounts(compressedDistanceSignal,numValues);
    
    probVector = AccSamp./(ones(numValues,1).*length(compressedDistanceSignal))';
    [dict,avglen] = huffmandictMod(intervalVector',probVector);
    [comp] = huffmanencoMod(compressedDistanceSignal,dict,intervalVector);
    signalSize = length(comp);
end