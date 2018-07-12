function [error,avglen,signalSize] = HuffmanDynamicPolar(signal,numValuesAng,numValuesRad,trueValue,plots,huffman)

avglen=ceil(log2(numValuesAng * numValuesRad));
signalSize = avglen*length(signal);

maxR = .65;
maxA = 2*pi;
startSig = zeros(length(signal),1);
distanceTotal = zeros(length(signal)-1,2);
startSig(1) = signal(1);

TreuValuesVecIndx = 1:trueValue:length(signal)-trueValue;

intervalRad = maxR/(numValuesRad-1);
intervalRadVector = [0:intervalRad:maxR]';
intervalRadVector2 = [0 + intervalRad/2:intervalRad:maxR - intervalRad/2]';
VectorIntervalsRad = intervalVariable(intervalRadVector2);
VectorIntervalsRad(end,2) = maxR - VectorIntervalsRad(end,1);

intervalAng = maxA/(numValuesAng-1);
intervalAngVector = [0:intervalAng:maxA]';
intervalAngVector2 = [intervalAng/2:intervalAng:maxA - intervalAng/2]';
VectorIntervalsAng = intervalVariable(intervalAngVector2);
VectorIntervalsAng(end,2) = maxA - VectorIntervalsAng(end,1);

intervalVector = intervalVectorFun(VectorIntervalsRad,VectorIntervalsAng);
intervalVector = round(intervalVector,6);

for i=1:length(signal)-1
    
    if(sum(TreuValuesVecIndx == (i+1)))
        startSig(i+1) = signal(i+1);
    else
        diffR = abs(signal(i+1) - startSig(i));
        diffA = atan2(imag(signal(i+1) - startSig(i)),real(signal(i+1) - startSig(i)));
        diffA = 2*pi*(diffA < 0) + diffA; 
        [~,indxR] = min(abs(diffR-intervalRadVector));
        [~,indxA] = min(abs(diffA-intervalAngVector));
        SAng = intervalAngVector(indxA);
        SRad = intervalRadVector(indxR);
        SAng = SAng - 2*pi*(SAng >= 2*pi);
        startSig(i+1) = startSig(i) + (SRad * cos(SAng) + 1i * SRad * sin(SAng));
        distanceTotal(i,:) = round([intervalRadVector(indxR),intervalAngVector(indxA)],6);
    end
    
end

distancePolar = abs(startSig(1:end-1))-abs(startSig(2:end));
distanceAngle = atan2(imag(startSig(1:end-1)),real(startSig(1:end-1))) - atan2(imag(startSig(2:end)),real(startSig(2:end)));


tmwaveformC = round(startSig,6);

error = EVM(signal,tmwaveformC,plots);

if(plots)
    
    figure
    histogram(distancePolar,intervalRadVector)
    xlabel('Intervals')
    ylabel('Num. samples')
    title('Sample difference distribution(Radius)')
    
    figure
    histogram(distanceAngle,intervalAngVector)
    xlabel('Intervals [radiants]')
    ylabel('Num. samples')
    title('Sample difference distribution(Angles)')
    
    figure
    plot(intervalVector(:,1),intervalVector(:,2), 'xr')
    title('Constellation')
    ylabel('Angles(rad)')
    xlabel('Radius')
    figure
    voronoi(intervalVector(:,1),intervalVector(:,2))
    title('Constellation')
    ylabel('Angles(rad)')
    xlabel('Radius')
    grid on
    
    figure
    plot(tmwaveformC, 'x')
    title('Signal compressed')
    ylabel('Quadrature')
    xlabel('Phase')
    
    intervalVectorTog = intervalVector(:,1) + 1i * intervalVector(:,2);
    
    vectorSectionSize = ones(numValuesRad,1)*(numValuesAng);    
    
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

    modulatedSignal = signalModulation(distanceTotal,intervalVector);
    AccSamp = contOcurrence(1:length(intervalVector),modulatedSignal);
    
    probVector = AccSamp./(ones(length(intervalVector),1).*length(distanceTotal))';
    [dict,avglen] = huffmandictMod(0:length(intervalVector)-1,probVector);
    [comp] = huffmanencoMod(modulatedSignal,dict,1:length(intervalVector));
    signalSize = length(comp);
end