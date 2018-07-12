function [error,avglen,signalSize] = HuffmanFitSplit(signal,numValues,plots,huffman)
avglen = ceil(log2(numValues));
signalSize = avglen*length(signal)*2;

numSamplesInt = ceil(2*length(signal)/numValues);
tmw = [real(signal);imag(signal)];
distancesSamples = sort(real(tmw));
pointsImp = distancesSamples(1:numSamplesInt:end);
pointsImp = [-1;pointsImp(2:end);1];
pointsImp = (pointsImp(1:end-1)+pointsImp(2:end))/2;

SignalPhase = real(signal);
SignalQuadrature = imag(signal);

constellation = pointsImp + 1i * pointsImp';
constellation = reshape(constellation,[],1);

minDAll = abs(signal-constellation');
[~,minInd] = min(minDAll,[],2);
tmwaveformC = constellation(minInd)'.';
error = EVM(signal,tmwaveformC,plots);
%% Plots
if(plots)
    figure
    plot(real(constellation),imag(constellation), 'xr')
    xlabel('Phase')
    ylabel('Quadrature')
    title('Constellation')
    figure
    voronoi(real(constellation),imag(constellation))
    title('Constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    
    figure
    plot(pointsImp,'x')
    title('Levels Assigments')
    ylabel('Interval')
    
    figure
    histogram(SignalPhase,pointsImp')
    title('Samples acumulated on each level (Phase)')
    ylabel('Samples accumulated')
    xlabel('Interval')
    
    figure
    histogram(SignalQuadrature,pointsImp')
    title('Samples acumulated on each level (Quadrature)')
    ylabel('Samples accumulated')
    xlabel('Interval')
    
    figure
    plot(tmwaveformC, 'x')
    title('Signal compressed')
    ylabel('Quadrature')
    xlabel('Phase')
    
    vectorSectionSize = ones(numValues,1)*numValues;
    
    errorsVector = [abs(tmwaveformC - signal),minInd];
    errorsHist = zeros(length(constellation),1);
    for i=1:length(constellation)
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
    
    tmwT = [real(signal);imag(signal)];
    minDAll = abs(tmwT-pointsImp');
    [~,minInd] = min(minDAll,[],2);
    tmwTotal = pointsImp(minInd)'.';
    
    histoEdges = 0.5:numValues+0.5;
    AccSamp = histcounts(minInd,histoEdges);
    
    probVector = AccSamp./(ones(numValues,1).*2*length(signal))';
    [dict,avglen] = huffmandictMod(round(pointsImp,5)',probVector);
    comp = huffmanencoMod(round(tmwTotal,5),dict,round(pointsImp,5));
    signalSize = length(comp);

end