function [error,avglen,signalSize] = HuffmanFitIQTogether(signal,numPhase,numQuadrature,plots,huffman)
avglen = 0;
signalSize = 0;

numSamplesInt = ceil(length(signal)/(numPhase));
distancesSamples = sort(real(signal));
pointsImpPhase = distancesSamples(1:numSamplesInt:end);
pointsImpPhase = [-1;pointsImpPhase(2:end);1];
pointsImpPhase = (pointsImpPhase(1:end-1)+pointsImpPhase(2:end))/2;

numSamplesInt = ceil(length(signal)/(numQuadrature));
distancesSamples = sort(imag(signal));
pointsImpQuad = distancesSamples(1:numSamplesInt:end);
pointsImpQuad = [-1;pointsImpQuad(2:end);1];
pointsImpQuad = (pointsImpQuad(1:end-1)+pointsImpQuad(2:end))/2;

constellation = pointsImpPhase + 1i * pointsImpQuad';
constellation = reshape(constellation,[],1);

minDAll = abs(signal-constellation');
[~,minInd] = min(minDAll,[],2);
tmwavesformC = constellation(minInd)'.';

error = EVM(signal,tmwavesformC,plots);

%% Plots with extra information about the compression
if(plots)
    
    voronoi(real(constellation),imag(constellation))
    title('Constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    axis([-1.1 1.1 -1.1 1.1])
    figure
    plot(real(constellation),imag(constellation), 'xr')
    xlabel('Phase')
    ylabel('Quadrature')
    title('Constellation')
    
    figure
    histogram(real(signal),pointsImpPhase')
    title('Samples acumulated on each level (Phase)')
    ylabel('Samples accumulated')
    xlabel('Interval')
    
    figure
    histogram(imag(signal),pointsImpQuad')
    title('Samples acumulated on each level (Quadrature)')
    ylabel('Samples accumulated')
    xlabel('Interval')
    
    figure
    plot(tmwavesformC, 'x')
    title('Signal compressed')
    ylabel('Quadrature')
    xlabel('Phase')
    
    histoEdges = 0.5:length(constellation)+0.5;
    
    vectorSectionSize = ones(numQuadrature,1)*numPhase;
    AccSamp = histcounts(minInd,histoEdges);
    modulatedSignalCell = mat2cell(AccSamp,1,vectorSectionSize');
    modulatedSignalTotal = cellfun(@(x) sum(x),modulatedSignalCell);
    
    figure
    bar(modulatedSignalTotal)
    title('Samples acumulated on each level')
    xlabel('Interval')
    ylabel('Samples acumulated')

    errorsVector = [abs(tmwavesformC - signal),minInd];
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

       
    histoEdges = 0.5:length(constellation)+0.5;
    AccSamp = histcounts(minInd,histoEdges);
    
    probVector = AccSamp./(ones(length(constellation),1).*length(signal))';
    [dict,avglen] = huffmandictMod(1:length(constellation),probVector);
    comp = huffmanencoMod(minInd,dict,1:length(constellation));
    signalSize = length(comp);

end