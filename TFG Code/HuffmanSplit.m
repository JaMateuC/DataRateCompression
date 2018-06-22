function [error,avglen,signalSize] = HuffmanSplit(signal,numBits,plots,huffman)
avglen = 0;
signalSize = 0;
maxI = 1;
 intervalAmplitude = 2*maxI/(numBits-1);
%   x = 0.5:numBits/2;
%   intervalAmplitudeVector = maxI* x.^1.45/max(x.^1.45);
%   intervalAmplitudeVector = [-fliplr(intervalAmplitudeVector),intervalAmplitudeVector]';
 
 intervalAmplitudeVector = [-maxI + intervalAmplitude/2:intervalAmplitude:...
    maxI-intervalAmplitude/2]';
VectorIntervals = intervalVariable(intervalAmplitudeVector);
VectorIntervals(end,2) = maxI - VectorIntervals(end,1);

tmwaveform2 = normalization(signal);

SignalPhase = real(tmwaveform2);
SignalQuadrature = imag(tmwaveform2);

compressedSignalPhase = signalCompression2(SignalPhase,VectorIntervals,maxI,-maxI);
compressedSignalQuadrature = signalCompression2(SignalQuadrature,...
    VectorIntervals,maxI,-maxI);
compressedSignal = round([compressedSignalPhase;compressedSignalQuadrature],5);
tmwaveformC = compressedSignalPhase + 1i * compressedSignalQuadrature;

error = EVM(tmwaveform2,tmwaveformC,plots);

%% Plots

if(plots)
    
    intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
    figure
    plot(intervalVector(:,1),intervalVector(:,2), 'xr')
    xlabel('Phase')
    ylabel('Quadrature')
    title('Constellation')
    
    plotSignal(signal)
    
    figure
    plot(intervalAmplitudeVector,'x')
    title('Levels Assigments')
    ylabel('Interval')
    
    figure
    histogram(SignalPhase,intervalAmplitudeVector')
    title('Samples acumulated on each level (Phase)')
    ylabel('Samples accumulated')
    xlabel('Interval')
    
    figure
    histogram(SignalQuadrature,intervalAmplitudeVector')
    title('Samples acumulated on each level (Quadrature)')
    ylabel('Samples accumulated')
    xlabel('Interval')
    
    figure
    plot(tmwaveformC, 'x')
    title('Signal compressed')
    ylabel('Quadrature')
    xlabel('Phase')
    
    ErrorAccumulated(compressedSignalPhase,compressedSignalQuadrature...
        ,intervalVector,tmwaveformC,tmwaveform2,VectorIntervals,VectorIntervals)
    
end

if(huffman)
    
    intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
    intervalVector = round(intervalVector(1:numBits,2),5);
    AccSamp = histcounts(compressedSignal,numBits);
    
    probVector = AccSamp./(ones(numBits,1).*2*length(signal))';
    [dict,avglen] = huffmandictMod(intervalVector',probVector);
    comp = huffmanencoMod(compressedSignal,dict,intervalVector);
    signalSize = length(comp);

end