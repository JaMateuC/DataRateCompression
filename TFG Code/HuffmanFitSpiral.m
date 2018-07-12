function [error,avglen,signalSize,newLen] = HuffmanFitSpiral(signal,numValues,numLoops,plots,huffman)
avglen = ceil(log2(numValues));
signalSize = avglen*length(signal);
maxI = numLoops*2*pi;
sens = 10000;

intervalAngle = maxI/(sens-1);
intervalAngle = 0:intervalAngle:maxI;
intervalRadius = 1/(numLoops*2*pi) * intervalAngle;
 
pointsQuant = intervalRadius.*cos(intervalAngle) + 1i * intervalRadius.*sin(intervalAngle);

minDAll = abs(signal-pointsQuant);
[~,minInd] = min(minDAll,[],2);

numSamplesInt = round(length(signal)/numValues);
AccSamp = histcounts(minInd,sens);
interv = 1;
j = 1;
z = 2;
while(z <= numValues && j <= sens)
    count = 0;
    while(count < numSamplesInt && j <= sens)
        count = count + AccSamp(j);
        j = j + 1;
    end
    interv(z) = j - 1;
    z = z + 1;
end

newQuant = pointsQuant(interv);
newQuant = (newQuant(1:end-1)+newQuant(2:end))/2;

minDAll = abs(signal-newQuant);
[~,minInd] = min(minDAll,[],2);
tmwaveformCFit = newQuant(minInd).';
AccSamp = histcounts(minInd,numValues);
newLen = length(newQuant);

error = EVM(signal,tmwaveformCFit,plots);

if(plots)
    
    figure
    plot(newQuant,'xr');
    title('Compressed signal constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    axis([-1 1 -1 1])
    
    figure
    voronoi(real(newQuant),imag(newQuant));
    title('Compressed signal constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    axis([-1 1 -1 1])
    
    figure
    bar(AccSamp)
    title('Samples acumulated on each value')
    xlabel('Interval')
    ylabel('Samples acumulated')
    
    figure
    plot(tmwaveformCFit, 'x')
    title('I/Q signal after compression')
    xlabel('Phase')
    ylabel('Quadrature')
    axis([-1 1 -1 1])
    
    vectorSectionSize = ones(length(newQuant),1);
    errorsVector = [abs(tmwaveformCFit - signal),minInd];
    errorsHist = zeros(length(newQuant),1);

    for i=1:length(newQuant)
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
    
    histoEdges = 0.5:newLen+0.5;
    AccSamp = histcounts(minInd,histoEdges);
    
    probVector = AccSamp./(ones(newLen,1).*length(signal)).';
    [dict,avglen] = huffmandictMod(round(newQuant,5),probVector);
    comp = huffmanencoMod(round(tmwaveformCFit,5),dict,round(newQuant,5).');
    signalSize = length(comp);

end