function [error] = HuffmanPolar(signal,numRadius,numAngles,plots)

% polar signal
tmwaveform2 = normalization(signal);
polarSignal = toPolar(tmwaveform2);

intervalAngle = 2*pi/numAngles;
intervalAngleVector = [0 + intervalAngle/2:intervalAngle:2*pi-intervalAngle/2]';
intervalRadius = 1/numRadius;
intervalRadiusVector = [0 + intervalRadius/2:intervalRadius:1-intervalRadius/2]';
intervalVectorRadius = intervalVariable(intervalRadiusVector);
intervalVectorRadius(end,2) = 1 - intervalVectorRadius(end,1);
intervalVectorAngle = intervalVariable(intervalAngleVector);
intervalVectorAngle(end,2) = 2*pi - intervalVectorAngle(end,1);

compressedSignal(:,1) = signalCompression(polarSignal(:,1),intervalVectorRadius,1,0);
compressedSignal(:,2) = signalCompression(polarSignal(:,2),intervalVectorAngle,0,0);
compressedSignal(:,2) = compressedSignal(:,2) .* ~(compressedSignal(:,2) >= 2*pi);

tmwavesformP = compressedSignal(:,1).*cos(compressedSignal(:,2)) + 1i * compressedSignal(:,1).*sin(compressedSignal(:,2));

error = EVM(tmwaveform2,tmwavesformP,plots);

%% Plots with extra information about the compression
if(plots)
    intervalVector = intervalVectorFun(intervalVectorRadius,intervalVectorAngle);
    intervalVectorInv = intervalVectorFun(intervalVectorAngle,intervalVectorRadius);

    constellation = intervalVector(:,1).*cos(intervalVector(:,2)) + 1i * intervalVector(:,1).*sin(intervalVector(:,2));
    figure
    plot(constellation,'xr');
    title(['Compressed signal constellation: Angle interval(' num2str(intervalRadius) '), Radius interval (' num2str(intervalAngle) ')'])
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
        ,intervalVector,tmwavesformP,tmwaveform2,intervalVectorRadius,intervalVectorAngle)
    ErrorAccumulated(compressedSignal(:,2),compressedSignal(:,1)...
        ,intervalVectorInv,tmwavesformP,tmwaveform2,intervalVectorAngle,intervalVectorRadius)

    figure
    tmwavesformP = compressedSignal(:,1).*cos(compressedSignal(:,2)) + 1i * compressedSignal(:,1).*sin(compressedSignal(:,2));
    plot(tmwavesformP, 'x')
    title('I/Q signal after compression')
    xlabel('Phase')
    ylabel('Quadrature')
end