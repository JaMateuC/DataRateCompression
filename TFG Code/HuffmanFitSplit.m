function [error,avglen,signalSize] = HuffmanFitSplit(signal,numBits,plots,huffman)
avglen = 0;
signalSize = 0;

tmwaveform2 = normalization(signal);

numSamplesInt = ceil(2*length(tmwaveform2)/numBits);
tmw = [real(tmwaveform2);imag(tmwaveform2)];
distancesSamples = sort(real(tmw));
pointsImp = distancesSamples(1:numSamplesInt:end);
pointsImp = [-1;pointsImp(2:end);1];
pointsImp = (pointsImp(1:end-1)+pointsImp(2:end))/2;

SignalPhase = real(tmwaveform2);
SignalQuadrature = imag(tmwaveform2);

constellation = pointsImp + 1i * pointsImp';
constellation = reshape(constellation,[],1);

minDAll = abs(tmwaveform2-constellation');
[~,minInd] = min(minDAll,[],2);
tmwaveformC = constellation(minInd)'.';
error = EVM(tmwaveform2,tmwaveformC,plots);
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
    
    vectorSectionSize = ones(numBits,1)*numBits;
    
    errorsVector = [abs(tmwaveformC - tmwaveform2),minInd];
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
    
    tmwT = [real(tmwaveform2);imag(tmwaveform2)];
    minDAll = abs(tmwT-pointsImp');
    [~,minInd] = min(minDAll,[],2);
    tmwTotal = pointsImp(minInd)'.';
    
    histoEdges = 0.5:numBits+0.5;
    AccSamp = histcounts(minInd,histoEdges);
    
    probVector = AccSamp./(ones(numBits,1).*2*length(signal))';
    [dict,avglen] = huffmandictMod(round(pointsImp,5)',probVector);
    comp = huffmanencoMod(round(tmwTotal,5),dict,round(pointsImp,5));
    signalSize = length(comp);

end