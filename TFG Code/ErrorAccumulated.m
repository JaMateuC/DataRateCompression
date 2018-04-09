function errorsHistTotal = ErrorAccumulated(compressedSignal1,compressedSignal2,intervalVector,tmwaveformCompressed,tmwaveformNormalized,VectorIntervals,VectorIntervals2)

    modulatedSignal = signalModulation([compressedSignal1,compressedSignal2],intervalVector); 
    errorsVector = [abs(tmwaveformCompressed - tmwaveformNormalized),modulatedSignal];
    errorsHist = zeros(length(intervalVector),1);
    vectorSectionSize = ones(length(VectorIntervals)+1,1)*(length(VectorIntervals2)+1);

    for i=1:length(intervalVector)
        errorsHist(i) = sum(errorsVector(errorsVector(:,2) == i));
    end
    errorsHistCell = mat2cell(errorsHist',1,vectorSectionSize');
    errorsHistTotal = cellfun(@(x) sum(x)/sum(errorsHist)*100,errorsHistCell);
    figure
    bar(errorsHistTotal)
    title('Error acumulation on each level')
    ylabel('Error acumulated')
    xlabel('Interval')