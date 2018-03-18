function [error] = HuffmanPolar(signal,numRadius,numAngles,plots)

% polar signal
tmwaveform2 = normalization(signal);
polarSignal = toPolar(tmwaveform2);

intervalAngle = 2*pi/numAngles;
intervalRadius = 1/numRadius;
intervalVectorRadius = intervalVariable(intervalRadius,0,1);
intervalVectorAngle = intervalVariable(intervalAngle,0,2*pi);

compressedSignal(:,1) = signalCompression(polarSignal(:,1),intervalVectorRadius,intervalRadius);
compressedSignal(:,2) = signalCompression(polarSignal(:,2),intervalVectorAngle,intervalAngle);
compressedSignal(:,2) = compressedSignal(:,2) .* ~(compressedSignal(:,2) >= 2*pi);

tmwavesformP = compressedSignal(:,1).*cos(compressedSignal(:,2)) + 1i * compressedSignal(:,1).*sin(compressedSignal(:,2));

error = EVM(tmwaveform2,tmwavesformP);

%% Inteting plots with extra information about the compression
if(plots)
    intervalVector = intervalVectorFun(intervalVectorRadius,intervalVectorAngle,intervalRadius,intervalAngle);
    intervalVectorInv = intervalVectorFun(intervalVectorAngle,intervalVectorRadius,intervalAngle,intervalRadius);

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
    figure
    histogram(modulatedSignalInv,histoEdges)


    histoEdges = 0.5:length(intervalVector)+0.5;
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
    hold off
    figure
    tmwavesformP = compressedSignal(:,1).*cos(compressedSignal(:,2)) + 1i * compressedSignal(:,1).*sin(compressedSignal(:,2));
    plot(tmwavesformP, 'x')
end