function [compressedSignal] = signalCompression2(signalInput,intervalVector,maxV,minV)

compressedSignal = zeros(length(signalInput),1);

for j=1:size(intervalVector,1)-1
    compressedSignal = (intervalVector(j,1) + intervalVector(j,2)).*(signalInput > intervalVector(j,1) & signalInput <= intervalVector(j+1,1)) +...
        ~(signalInput > intervalVector(j,1) & signalInput <= intervalVector(j+1,1)).*compressedSignal;
end

compressedSignal = maxV.*(signalInput > intervalVector(end,1)) + ~(signalInput > intervalVector(end,1)).*compressedSignal;
compressedSignal = minV.*(signalInput <= intervalVector(1,1)) + ~(signalInput <= intervalVector(1,1)).*compressedSignal;