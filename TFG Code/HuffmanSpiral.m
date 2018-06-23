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
    plot(minDAll(minInd))
    title('Error of each sample')
    xlabel('Sample')
    ylabel('Euclidean error distance')
    
    figure
    plot(tmwaveformC, 'x')
    title('I/Q signal after compression')
    xlabel('Phase')
    ylabel('Quadrature')
    axis([-1 1 -1 1])
end
if(huffman)
    
    
    AccSamp = histcounts(minInd,numBits);
    
    probVector = AccSamp./(ones(numBits,1).*length(signal)).';
    [dict,avglen] = huffmandictMod(pointsQuant,probVector);
    comp = huffmanencoMod(tmwaveformC,dict,pointsQuant.');
    signalSize = length(comp);

end