function [error,avglen,signalSize] = HuffmanSplit(signal,numValues,plots,huffman)
avglen = ceil(log2(numValues));
signalSize = avglen*2*length(signal);
maxI = 1;
 intervalAmplitude = 2*maxI/(numValues-1);
 
 intervalAmplitudeVector = [-maxI + intervalAmplitude/2:intervalAmplitude:...
    maxI-intervalAmplitude/2]';
VectorIntervals = intervalVariable(intervalAmplitudeVector);
VectorIntervals(end,2) = maxI - VectorIntervals(end,1);

SignalPhase = real(signal);
SignalQuadrature = imag(signal);

compressedSignalPhase = signalCompression2(SignalPhase,VectorIntervals,maxI,-maxI);
compressedSignalQuadrature = signalCompression2(SignalQuadrature,...
    VectorIntervals,maxI,-maxI);
if(sum(real(signal)) == 0)
    compressedSignalPhase = compressedSignalPhase .* 0;
end
if(sum(imag(signal)) == 0)
    compressedSignalQuadrature = compressedSignalQuadrature .* 0;
end
compressedSignal = round([compressedSignalPhase;compressedSignalQuadrature],6);
tmwaveformC = compressedSignalPhase + 1i * compressedSignalQuadrature;

error = EVM(signal,tmwaveformC,plots);

%% Plots

if(plots)
    
    intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
    figure
    plot(intervalVector(:,1),intervalVector(:,2), 'xr')
    xlabel('Phase')
    ylabel('Quadrature')
    title('Constellation')
    
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
        ,intervalVector,tmwaveformC,signal,VectorIntervals,VectorIntervals)
    
end

if(huffman)
    intervalVector = intervalVectorFun(VectorIntervals,VectorIntervals);
    intervalVector = round(intervalVector(1:numValues,2),6);
    
    if(sum(imag(signal)) == 0)
        AccSamp = contOcurrence(intervalVector,compressedSignal(1:length(signal)));
        probVector = AccSamp./(ones(numValues,1).*length(signal))';
        [dict,avglen] = huffmandictMod(intervalVector',probVector);
        comp = huffmanencoMod(compressedSignal(1:length(signal)),dict,intervalVector);
        signalSize = length(comp);
    elseif(sum(real(signal)) == 0)
        AccSamp = contOcurrence(intervalVector,compressedSignal(length(signal):end));
        probVector = AccSamp./(ones(numValues,1).*length(signal))';
        [dict,avglen] = huffmandictMod(intervalVector',probVector);
        comp = huffmanencoMod(compressedSignal(length(signal):end),dict,intervalVector);
        signalSize = length(comp);
    else

        AccSamp = contOcurrence(intervalVector,compressedSignal);
        probVector = AccSamp./(ones(numValues,1).*2*length(signal))';
        [dict,avglen] = huffmandictMod(intervalVector',probVector);
        comp = huffmanencoMod(compressedSignal,dict,intervalVector);
        signalSize = length(comp);
    end  

end