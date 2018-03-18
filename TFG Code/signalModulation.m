function [modulatedSignal] = signalModulation(compressedSignal,intervalVector)

modulatedSignal = zeros(length(compressedSignal),1);

for i=1:length(modulatedSignal)
    
    for j=1:length(intervalVector)
        if(round(compressedSignal(i,:),7) == round(intervalVector(j,:),7))
            modulatedSignal(i) = j;
            break
        end
    end
    
end