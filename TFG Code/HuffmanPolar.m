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

compressedSignal(:,1) = signalCompression(polarSignal(:,1),intervalVectorRadius,1);
compressedSignal(:,2) = signalCompression(polarSignal(:,2),intervalVectorAngle,0);
compressedSignal(:,2) = compressedSignal(:,2) .* ~(compressedSignal(:,2) >= 2*pi);

tmwavesformP = compressedSignal(:,1).*cos(compressedSignal(:,2)) + 1i * compressedSignal(:,1).*sin(compressedSignal(:,2));

error = EVM(tmwaveform2,tmwavesformP);

%% Inteting plots with extra information about the compression
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

    modulatedSignal = signalModulation(compressedSignal,intervalVector);

    compressedSignalInv = [compressedSignal(:,2),compressedSignal(:,1)];
    modulatedSignalInv = signalModulation(compressedSignalInv,intervalVectorInv);

    histoEdges = 0.5:length(intervalVector)+0.5;

    figure
    histogram(modulatedSignal,histoEdges)
    title('Samples acumulated on each level (Radius intervals)')
    xlabel('Levels')
    ylabel('Samples acumulated')
    figure
    histogram(modulatedSignalInv,histoEdges)
    title('Samples acumulated on each level (Angles intervals)')
    xlabel('Levels')
    ylabel('Samples acumulated')

    nCounts = histcounts(modulatedSignal,histoEdges);
    totalCount = sum(nCounts);
    modulatedProbability = nCounts/totalCount;
    dict = huffmandict(1:length(intervalVector),modulatedProbability);

    errorsVector = [abs(tmwavesformP - tmwaveform2),modulatedSignal];
    errorsHist = zeros(length(dict),1);

    for i=1:length(dict)
        errorsHist(i) = sum(errorsVector(errorsVector(:,2) == i));
    end
    figure
    bar(errorsHist)
    hold on 
    x = 0:length(dict);
    y = raylpdf(x,150);
    y = (max(errorsHist)-1)/max(y) * y;
    plot(x,y,'g');
    title('Error acumulation on each level')
    ylabel('Error acumulated')
    xlabel('Levels')
    hold off
    figure
    tmwavesformP = compressedSignal(:,1).*cos(compressedSignal(:,2)) + 1i * compressedSignal(:,1).*sin(compressedSignal(:,2));
    plot(tmwavesformP, 'x')
    title('I/Q signal after compression')
    xlabel('Phase')
    ylabel('Quadrature')
end