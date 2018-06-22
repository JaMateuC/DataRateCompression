function [error,avglen,signalSize] = HuffmanFitPolar(polarSignal,normSignal,numRadius,numAngles,plots,huffman)
avglen = 0;
signalSize = 0;

% polar signal
maxR = 2*pi;

numSamplesInt = ceil(length(polarSignal)/(numRadius));
distancesSamples = sort(abs(normSignal));
pointsImp = distancesSamples(1:numSamplesInt:end);
pointsImp = [0;pointsImp(2:end);1];

pointsImp = (pointsImp(1:end-1)+pointsImp(2:end))/2;

intervalAngle = maxR/numAngles;
intervalAngleVector = [0 + intervalAngle/2:intervalAngle:maxR-intervalAngle/2]';

intervalVector = pointsImp + 1i * intervalAngleVector';
intervalVector = reshape(intervalVector,[],1);
intervalVector2(:,2) = imag(intervalVector);
intervalVector2(:,1) = real(intervalVector);
constellation = intervalVector2(:,1).*cos(intervalVector2(:,2)) + 1i * intervalVector2(:,1).*sin(intervalVector2(:,2));
constellation = real(constellation) + 1i * round(imag(constellation),5);
constellation = unique(constellation);

minDAll = abs(normSignal-constellation');
[~,minInd] = min(minDAll,[],2);
tmwavesformP = constellation(minInd)'.';

error = EVM(normSignal,tmwavesformP,plots);

%% Plots with extra information about the compression
if(plots)

    figure
    voronoi(real(constellation),imag(constellation))
    title('Constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on
    figure
    plot(constellation,'xr');
    title('Constellation')
    xlabel('Phase')
    ylabel('Quadrature')
    grid on

    histoEdges = 0.5:length(intervalVector2)+0.5;
    
    vectorSectionSize = ones(numRadius,1)*length(intervalAngleVector);
%     vectorSectionSize = [1;vectorSectionSize];
    AccSamp = histcounts(minInd,histoEdges);
    modulatedSignalCell = mat2cell(AccSamp,1,vectorSectionSize');
    modulatedSignalTotal = cellfun(@(x) sum(x),modulatedSignalCell);
    
    figure
    bar(modulatedSignalTotal)
    title('Samples acumulated on each level (Radius intervals)')
    xlabel('Interval')
    ylabel('Samples acumulated')

    errorsVector = [abs(tmwavesformP - normSignal),minInd];
    errorsHist = zeros(length(intervalVector2),1);

    for i=1:length(intervalVector2)
        errorsHist(i) = sum(errorsVector(errorsVector(:,2) == i));
    end
    errorsHistCell = mat2cell(errorsHist',1,vectorSectionSize');
    errorsHistTotal = cellfun(@(x) sum(x)/sum(errorsHist)*100,errorsHistCell);
    figure
    bar(errorsHistTotal)
    title('Error acumulation on each level')
    ylabel('Error acumulated')
    xlabel('Interval')

    figure
    plot(tmwavesformP, 'x')
    title('I/Q signal after compression')
    xlabel('Phase')
    ylabel('Quadrature')
end

if(huffman)
    
    histoEdges = 0.5:length(intervalVector)+0.5;
    AccSamp = histcounts(minInd,histoEdges);
    
    probVector = AccSamp./(ones(length(intervalVector),1).*length(polarSignal(:,1)))';
    [dict,avglen] = huffmandictMod(1:length(intervalVector),probVector);
    comp = huffmanencoMod(minInd,dict,1:length(intervalVector));
    signalSize = length(comp);

end