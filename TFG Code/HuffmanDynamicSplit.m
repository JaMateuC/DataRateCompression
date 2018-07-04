function [error,avglen,signalSize] = HuffmanDynamicSplit(signal,numValues,trueValue,plots,huffman)

avglen = 0;
signalSize = 0;

maxI = .7;
startSigP = zeros(length(signal),1);
startSigQ = zeros(length(signal),1);

startSigP(1) = real(signal(1));
startSigQ(1) = imag(signal(1));

distanceReal = real(signal(1:end-1))-real(signal(2:end));
distanceImag = imag(signal(1:end-1))-imag(signal(2:end));

TreuValuesVecIndx = 1:trueValue:length(signal)-trueValue;

intervalAmplitude = 2*maxI/(numValues-1);
intervalAmplitudeVector = [-maxI + intervalAmplitude/2:intervalAmplitude:...
    maxI-intervalAmplitude/2]';
VectorIntervals = intervalVariable(intervalAmplitudeVector);
VectorIntervals(end,2) = maxI - VectorIntervals(end,1);
intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
const = intervalVector(1:numValues:end,1);

j = 1;

for i=1:length(distanceImag)
    
    if(sum(TreuValuesVecIndx == (i+1)))
        startSigP(i+1) = real(signal(i+1));
        startSigQ(i+1) = imag(signal(i+1));
    else
        diffP = real(signal(i+1)) - startSigP(i);
        diffQ = imag(signal(i+1)) - startSigQ(i);
        [~,indxP] = min(abs(diffP-const));
        [~,indxQ] = min(abs(diffQ-const));
        compressedDistanceSignalPhase(j) = const(indxP);
        compressedDistanceSignalQuadrature(j) = const(indxQ);
        startSigP(i+1) = startSigP(i) + const(indxP);
        startSigQ(i+1) = startSigQ(i) + const(indxQ);
        j = j + 1;
    end
    
end

compressedDistanceSignal = round([compressedDistanceSignalPhase';compressedDistanceSignalQuadrature'],6);
tmwaveformC = startSigP + 1i * startSigQ;

error = EVM(signal,tmwaveformC,plots);

if(plots)
    
    figure
    histogram(compressedDistanceSignalPhase,intervalAmplitudeVector')
    xlabel('Intervals')
    ylabel('Num. samples')
    title('Sample difference distribution(Phase)')
    
    figure
    histogram(compressedDistanceSignalQuadrature,intervalAmplitudeVector')
    xlabel('Intervals')
    ylabel('Num. samples')
    title('Sample difference distribution(Quadrature)')
    
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
    
    minDAll = abs(signal-intervalVectorTog');
    [~,minInd] = min(minDAll,[],2);
    
    errorsVector = [abs(tmwaveformC - signal),minInd];
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

    const = round(const,6);
    AccSamp = contOcurrence(const,compressedDistanceSignal);
    
    probVector = AccSamp./(ones(numValues,1).*length(compressedDistanceSignal))';
    [dict,avglen] = huffmandictMod(const',probVector);
    [comp] = huffmanencoMod(compressedDistanceSignal,dict,const);
    signalSize = length(comp);
end