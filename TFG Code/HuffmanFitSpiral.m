function [error,avglen,signalSize,newLen] = HuffmanFitSpiral(signal,numBits,numV,plots,huffman)
avglen = 0;
signalSize = 0;
maxI = numV*2*pi;
sens = 10000;

intervalAngle = maxI/(sens-1);
intervalAngle = 0:intervalAngle:maxI;
intervalRadius = 1/(numV*2*pi) * intervalAngle;
 
pointsQuant = intervalRadius.*cos(intervalAngle) + 1i * intervalRadius.*sin(intervalAngle);

minDAll = abs(signal-pointsQuant);
[~,minInd] = min(minDAll,[],2);

numSamplesInt = round(length(signal)/numBits);
AccSamp = histcounts(minInd,sens);
interv = 1;
j = 1;
z = 2;
while(z <= numBits && j <= sens)
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
AccSamp = histcounts(minInd,numBits);
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
    plot(minDAll(minInd))
    title('Error of each sample')
    xlabel('Sample')
    ylabel('Euclidean error distance')
    
    figure
    plot(tmwaveformCFit, 'x')
    title('I/Q signal after compression')
    xlabel('Phase')
    ylabel('Quadrature')
    axis([-1 1 -1 1])
end

if(huffman)
    
    AccSamp = histcounts(minInd,newLen);
    
    probVector = AccSamp./(ones(newLen,1).*length(signal)).';
    [dict,avglen] = huffmandictMod(round(newQuant,5),probVector);
    comp = huffmanencoMod(round(tmwaveformCFit,5),dict,round(newQuant,5).');
    signalSize = length(comp);

end