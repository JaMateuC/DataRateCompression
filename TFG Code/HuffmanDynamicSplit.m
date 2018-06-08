function [error,avglen,signalSize] = HuffmanDynamicSplit(signal,numValues,trueValue,plots,huffman)

maxI = .7;
startSigP = zeros(length(signal),1);
startSigQ = zeros(length(signal),1);

tmwaveform2 = normalization(signal);
startSigP(1) = real(tmwaveform2(1));
startSigQ(1) = imag(tmwaveform2(1));

distanceReal = real(tmwaveform2(1:end-1))-real(tmwaveform2(2:end));
distanceImag = imag(tmwaveform2(1:end-1))-imag(tmwaveform2(2:end));

TreuValuesVecIndx = 1:trueValue:length(tmwaveform2)-trueValue;
TreuValuesVec = tmwaveform2(TreuValuesVecIndx);


intervalAmplitude = 2*maxI/(numValues-1);
 intervalAmplitudeVector = [-maxI + intervalAmplitude/2:intervalAmplitude:...
    maxI-intervalAmplitude/2]';
VectorIntervals = intervalVariable(intervalAmplitudeVector);
VectorIntervals(end,2) = maxI - VectorIntervals(end,1);

compressedDistanceSignalPhase = signalCompression2(distanceReal,VectorIntervals,maxI,-maxI);
compressedDistanceSignalQuadrature = signalCompression2(distanceImag,...
    VectorIntervals,maxI,-maxI);

for i=1:length(distanceImag)
    
    if(sum(TreuValuesVecIndx == (i+1)))
        startSigP(i+1) = real(tmwaveform2(i+1));
        startSigQ(i+1) = imag(tmwaveform2(i+1));
    else
        startSigP(i+1) = startSigP(i) - compressedDistanceSignalPhase(i);
        startSigQ(i+1) = startSigQ(i) - compressedDistanceSignalQuadrature(i);
    end
    
end

compressedDistanceSignal = round([compressedDistanceSignalPhase;compressedDistanceSignalQuadrature],5);
tmwaveformC = startSigP + 1i * startSigQ;

error = EVM(tmwaveform2,tmwaveformC,plots);

if(plots)
    
    figure
    histogram(distanceReal,200)

    figure
    histogram(distanceImag,200)
    
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