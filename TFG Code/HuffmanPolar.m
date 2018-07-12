function [error,avglen,signalSize] = HuffmanPolar(polarSignal,normSignal,numRadius,numAngles,plots,huffman)
avglen = 0;
signalSize = avglen*length(signal);

% polar signal
maxR = 2*pi;
maxA = 1;

intervalAngle = maxR/numAngles;
intervalAngleVector = [0 + intervalAngle/2:intervalAngle:maxR-intervalAngle/2]';
intervalRadius = maxA/numRadius;
intervalRadiusVector = [0 + intervalRadius/2:intervalRadius:maxA-intervalRadius/2]';
intervalVectorRadius = intervalVariable(intervalRadiusVector);
intervalVectorRadius(end,2) = maxA - intervalVectorRadius(end,1);
intervalVectorAngle = intervalVariable(intervalAngleVector);
intervalVectorAngle(end,2) = maxR - intervalVectorAngle(end,1);

compressedSignal = signalCompression(polarSignal,intervalVectorRadius,intervalVectorAngle,maxA,0,0,0);
compressedSignal(:,2) = compressedSignal(:,2) .* ~(compressedSignal(:,2) >= maxR);

tmwavesformP = compressedSignal(:,1).*cos(compressedSignal(:,2)) + 1i * compressedSignal(:,1).*sin(compressedSignal(:,2));

error = EVM(normSignal,tmwavesformP,plots);

%% Plots with extra information about the compression
if(plots)
    intervalVector = intervalVectorFun(intervalVectorRadius,intervalVectorAngle);
    intervalVectorInv = intervalVectorFun(intervalVectorAngle,intervalVectorRadius);

    constellation = intervalVector(:,1).*cos(intervalVector(:,2)) + 1i * intervalVector(:,1).*sin(intervalVector(:,2));
    constellation = real(constellation) + 1i * round(imag(constellation),5);
    constellation = unique(constellation);
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

    histoEdges = 0.5:length(intervalVector)+0.5;
    
    modulatedSignal = signalModulation(compressedSignal,intervalVector);
    vectorSectionSize = ones(length(intervalRadiusVector)+1,1)*(length(intervalAngleVector)+1);
    AccSamp = histcounts(modulatedSignal,histoEdges);
    modulatedSignalCell = mat2cell(AccSamp,1,vectorSectionSize');
    modulatedSignalTotal = cellfun(@(x) sum(x),modulatedSignalCell);
    
    figure
    bar(modulatedSignalTotal)
    title('Samples acumulated on each level (Radius intervals)')
    xlabel('Interval')
    ylabel('Samples acumulated')

    compressedSignalInv = [compressedSignal(:,2),compressedSignal(:,1)];
    modulatedSignalInv = signalModulation(compressedSignalInv,intervalVectorInv);
    vectorSectionSizeInv = ones(length(intervalAngleVector)+1,1)*(length(intervalRadiusVector)+1);
    AccSampInv = histcounts(modulatedSignalInv,histoEdges);
    modulatedSignalCellInv = mat2cell(AccSampInv,1,vectorSectionSizeInv');
    modulatedSignalTotalInv = cellfun(@(x) sum(x),modulatedSignalCellInv);

    figure
    bar(modulatedSignalTotalInv)
    title('Samples acumulated on each level (Angles intervals)')
    xlabel('Interval')
    ylabel('Samples acumulated')

    ErrorAccumulated(compressedSignal(:,1),compressedSignal(:,2)...
        ,intervalVector,tmwavesformP,normSignal,intervalVectorRadius,intervalVectorAngle)
    ErrorAccumulated(compressedSignal(:,2),compressedSignal(:,1)...
        ,intervalVectorInv,tmwavesformP,normSignal,intervalVectorAngle,intervalVectorRadius)

    figure
    plot(tmwavesformP, 'x')
    title('I/Q signal after compression')
    xlabel('Phase')
    ylabel('Quadrature')
end

%% Huffman
if(huffman)
    compressedSignal(:,2) = compressedSignal(:,2) .* ~(compressedSignal(:,1) == 0);
    intervalVector = intervalVectorFun(intervalVectorRadius(2:end,:),intervalVectorAngle(1:end-1,:));
    intervalVector = [[0,0];intervalVector];
    
    histoEdges = 0.5:length(intervalVector)+0.5;
    modulatedSignal = signalModulation(compressedSignal,intervalVector);
    AccSamp = histcounts(modulatedSignal,histoEdges);
    
    probVector = AccSamp./(ones(length(intervalVector),1).*length(polarSignal(:,1)))';
    [dict,avglen] = huffmandictMod(1:length(intervalVector),probVector);
    comp = huffmanencoMod(modulatedSignal,dict,1:length(intervalVector));
    signalSize = length(comp);

end