function [AccSamples] = contOcurrence(levels,compressedSignal)

AccSamples = zeros(1,length(levels));

for i=1:length(levels)
    AccSamples(i) = sum(levels(i) == compressedSignal);
end