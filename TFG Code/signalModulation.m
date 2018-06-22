function [modulatedSignal] = signalModulation(compressedSignal,intervalVector)

modulatedSignal = zeros(length(compressedSignal),1);
compressedSignal = round(compressedSignal,7);
intervalVector = round(intervalVector,7);

for j=1:length(intervalVector)
    modulatedSignal = j.*((compressedSignal(:,1) == intervalVector(j,1)) & (compressedSignal(:,2) == intervalVector(j,2)))+...
    ~((compressedSignal(:,1) == intervalVector(j,1)) & (compressedSignal(:,2) == intervalVector(j,2))).*modulatedSignal;
end