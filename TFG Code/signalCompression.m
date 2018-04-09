function [compressedSignal] = signalCompression(signalInput,intervalVector,max,min)

compressedSignal = zeros(length(signalInput),1);

for i=1:length(signalInput)
    if(signalInput(i) > intervalVector(end,1))
        compressedSignal(i) = max;
    elseif(signalInput(i) <= intervalVector(1,1))
        compressedSignal(i) = min;
    else
        for j=1:length(intervalVector)-1
            if(signalInput(i) > intervalVector(j,1) && signalInput(i) <= intervalVector(j+1,1))
                compressedSignal(i) = intervalVector(j,1) + intervalVector(j,2);
                break;
            end
        end
    end
end