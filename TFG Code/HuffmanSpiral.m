function [error,avglen,signalSize] = HuffmanSpiral(signal,numBits,numV,plots,huffman)
avglen = 0;
signalSize = 0;
maxI = numV*2*pi;

intervalAngle = maxI/(numBits-1);
intervalAngle = 0:intervalAngle:maxI;
intervalRadius = 1/(numV*2*pi) * intervalAngle;
 
pointsQuant = intervalRadius.*cos(intervalAngle) + 1i * intervalRadius.*sin(intervalAngle);

minDAll = abs(signal-pointsQuant);
[~,minInd] = min(minDAll,[],2);
tmwaveformC = pointsQuant(minInd).';
error = EVM(signal,tmwaveformC,plots);
if(plots)
    
    figure
    plot(pointsQuant,'xr');
    title('Compressed signal constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    axis([-1 1 -1 1])
    
    figure
    voronoi(real(pointsQuant),imag(pointsQuant));
    title('Compressed signal constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    axis([-1 1 -1 1])
    
    AccSamp = histcounts(minInd,numBits);
    
    figure
    bar(AccSamp)
    title('Samples acumulated on each value')
    xlabel('Interval')
    ylabel('Samples acumulated')
    
    figure
    plot(tmwaveformC, 'x')
    title('I/Q signal after compression')
    xlabel('Phase')
    ylabel('Quadrature')
    axis([-1 1 -1 1])
    
    vectorSectionSize = ones(length(pointsQuant),1);
    errorsVector = [abs(tmwaveformC - signal),minInd];
    errorsHist = zeros(length(pointsQuant),1);

    for i=1:length(pointsQuant)
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
    
    histoEdges = 0.5:numBits+0.5;
    AccSamp = histcounts(minInd,histoEdges);
    
    probVector = AccSamp./(ones(numBits,1).*length(signal)).';
    [dict,avglen] = huffmandictMod(pointsQuant,probVector);
    comp = huffmanencoMod(tmwaveformC,dict,pointsQuant.');
    signalSize = length(comp);

end