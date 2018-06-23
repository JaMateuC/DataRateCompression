function [error,avglen,signalSize] = HuffmanDynamicPolar(signal,numValuesAng,numValuesRad,trueValue,plots,huffman)

avglen=0;
signalSize = 0;

maxR = .65;
maxA = 2*pi;
startSig = zeros(length(signal),1);
startSig(1) = signal(1);

TreuValuesVecIndx = 1:trueValue:length(signal)-trueValue;

intervalRad = 2*maxR/(numValuesRad-1);
intervalRadVector = [-maxR:intervalRad:maxR]';
VectorIntervalsRad = intervalVariable(intervalRadVector);
VectorIntervalsRad(end,2) = maxR - VectorIntervalsRad(end,1);

intervalAng = maxA/(numValuesAng-1);
intervalAngVector = [0:intervalAng:maxA]';
VectorIntervalsAng = intervalVariable(intervalAngVector);
VectorIntervalsAng(end,2) = maxA - VectorIntervalsAng(end,1);
intervalVector = intervalVectorFun(VectorIntervalsRad,VectorIntervalsAng);

for i=1:length(signal)-1
    
    if(sum(TreuValuesVecIndx == (i+1)))
        startSig(i+1) = signal(i+1);
    else
        diffR = abs(signal(i+1)) - abs(startSig(i));
        diffA = atan2(imag(startSig(i)),real(startSig(i))) - atan2(imag(signal(i+1)),real(signal(i+1)));
        diffA = 2*pi*(diffA < 0) + diffA; 
        [~,indxR] = min(abs(diffR-intervalRadVector));
        [~,indxA] = min(abs(diffA-intervalAngVector));
        SAng = atan2(imag(startSig(i)),real(startSig(i))) - intervalAngVector(indxA);
        SRad = abs(startSig(i)) + intervalRadVector(indxR);
        SAng = SAng - 2*pi*(SAng >= 2*pi);
        startSig(i+1) = SRad * cos(SAng) + 1i * SRad * sin(SAng);
    end
    
end

distancePolar = abs(startSig(1:end-1))-abs(startSig(2:end));
distanceAngle = atan2(imag(startSig(1:end-1)),real(startSig(1:end-1))) - atan2(imag(startSig(2:end)),real(startSig(2:end)));

tmwaveformC = round(startSig,5);

error = EVM(signal,tmwaveformC,plots);

if(plots)
    
    figure
    histogram(distancePolar,100)
    
    figure
    histogram(distanceAngle,100)
    
    figure
    plot(intervalVector(:,1),intervalVector(:,2), 'xr')
    xlabel('Phase')
    ylabel('Quadrature')
    title('Constellation')
    figure
    voronoi(intervalRadVector,intervalAngVector)
    title('Constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    
    figure
    plot(intervalAmplitudeVector,'x')
    title('Levels Assigment')
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
    intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
    intervalVector = round(intervalVector(1:numValues,2),5);
    AccSamp = histcounts(compressedDistanceSignal,numValues);
    
    probVector = AccSamp./(ones(numValues,1).*length(compressedDistanceSignal))';
    [dict,avglen] = huffmandictMod(intervalVector',probVector);
    [comp] = huffmanencoMod(compressedDistanceSignal,dict,intervalVector);
    signalSize = length(comp);
end