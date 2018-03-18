function [error] = HuffmanSplit(signal,numBits,plots)

intervalAmplitude = 2/numBits;

VectorIntervals = intervalVariable(intervalAmplitude,-1,1);

tmwaveform2 = normalization(signal);

SignalPhase = real(tmwaveform2);
SignalQuadrature = imag(tmwaveform2);

compressedSignalPhase = signalCompression(SignalPhase,VectorIntervals,intervalAmplitude);
compressedSignalQuadrature = signalCompression(SignalQuadrature,VectorIntervals,intervalAmplitude);
tmwaveformC = compressedSignalPhase + 1i * compressedSignalQuadrature;

error = EVM(tmwaveform2,tmwaveformC);

if(plots)
end