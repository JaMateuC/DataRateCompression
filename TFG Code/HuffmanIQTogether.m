function [error] = HuffmanIQTogether(signal,numPhase,numQuadrature,plots)

% polar signal
tmwaveform2 = normalization(signal);

intervalPhase = 2/numPhase;
intervalPhaseVector = [-1 + intervalPhase/2:intervalPhase:1-intervalPhase/2]';
intervalQuadrature = 2/numQuadrature;
intervalQuadratureVector = [-1 + intervalQuadrature/2:intervalQuadrature:1-intervalQuadrature/2]';
intervalVectorQuadrature = intervalVariable(intervalQuadratureVector);
intervalVectorQuadrature(end,2) = 1 - intervalVectorQuadrature(end,1);
intervalVectorPhase = intervalVariable(intervalPhaseVector);
intervalVectorPhase(end,2) = 1 - intervalVectorPhase(end,1);

compressedSignal(:,1) = signalCompression2(real(tmwaveform2),intervalVectorPhase,1,-1);
compressedSignal(:,2) = signalCompression2(imag(tmwaveform2),intervalVectorQuadrature,1,-1);

tmwavesformC = compressedSignal(:,1) + 1i * compressedSignal(:,2);

error = EVM(tmwaveform2,tmwavesformC,plots);

%% Plots with extra information about the compression
if(plots)
    intervalVector = intervalVectorFun(intervalVectorPhase,intervalVectorQuadrature);
    figure
    plot(intervalVector(:,1),intervalVector(:,2), 'xr')
    xlabel('Phase')
    ylabel('Quadrature')
    title('Constellation')
    
    plotSignal(signal)
    
    figure
    histogram(real(tmwaveform2),intervalPhaseVector')
    title('Samples acumulated on each level (Phase)')
    ylabel('Samples accumulated')
    xlabel('Interval')
    
    figure
    histogram(imag(tmwaveform2),intervalQuadratureVector')
    title('Samples acumulated on each level (Quadrature)')
    ylabel('Samples accumulated')
    xlabel('Interval')
    
    figure
    plot(tmwavesformC, 'x')
    title('Signal compressed')
    ylabel('Quadrature')
    xlabel('Phase')
    
    ErrorAccumulated(compressedSignal(:,1),compressedSignal(:,2)...
        ,intervalVector,tmwavesformC,tmwaveform2,intervalVectorQuadrature,intervalVectorPhase)
end