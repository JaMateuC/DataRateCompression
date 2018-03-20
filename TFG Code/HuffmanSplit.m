function [error] = HuffmanSplit(signal,numBits,plots)

% intervalAmplitude = 2/numBits;
 x = 0.5:numBits/2;
 intervalAmplitudeVector = 0.9* x.^1.45/max(x.^1.45);
 intervalAmplitudeVector = [-fliplr(intervalAmplitudeVector),intervalAmplitudeVector]';
 
% intervalAmplitudeVector = [-1 + intervalAmplitude/2:intervalAmplitude:1-intervalAmplitude/2]';
VectorIntervals = intervalVariable(intervalAmplitudeVector);
VectorIntervals(end,2) = 1 - VectorIntervals(end,1);

tmwaveform2 = normalization(signal);

SignalPhase = real(tmwaveform2);
SignalQuadrature = imag(tmwaveform2);

compressedSignalPhase = signalCompression(SignalPhase,VectorIntervals,1);
compressedSignalQuadrature = signalCompression(SignalQuadrature,VectorIntervals,1);
tmwaveformC = compressedSignalPhase + 1i * compressedSignalQuadrature;

error = EVM(tmwaveform2,tmwaveformC);

if(plots)
    
    intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
    figure
    plot(intervalVector(:,1),intervalVector(:,2), 'xr')
    xlabel('Phase')
    ylabel('Quadrature')
    title('Constellation')
    
    figure
    plotSignal(signal)
    
    figure
    plot(intervalAmplitudeVector,'x')
    title('Levels Assigments')
    ylabel('intervals')
    
    figure
    histogram(SignalPhase,intervalAmplitudeVector')
    title('Samples acumulated on each level (Phase)')
    ylabel('Samples accumulated')
    xlabel('Levels')
    
    figure
    histogram(SignalQuadrature,intervalAmplitudeVector')
    title('Samples acumulated on each level (Quadrature)')
    ylabel('Samples accumulated')
    xlabel('Levels')
    
    figure
    plot(tmwaveformC, 'x')
    title('Signal compressed')
    ylabel('Quadrature')
    xlabel('Phase')
    
end