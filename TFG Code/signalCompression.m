function [compressedSignal] = signalCompression(signalInput,intervalVector,interval)

compressedSignal = zeros(length(signalInput),1);

for i=1:length(signalInput)
    if(signalInput(i) > intervalVector(end,1))
        compressedSignal(i) = intervalVector(end,1) + interval/2;
    elseif(signalInput(i) <= intervalVector(1,1))
        compressedSignal(i) = 0;
    else
        for j=1:length(intervalVector)-1
            if(signalInput(i) > intervalVector(j,1) && signalInput(i) <= intervalVector(j+1,1))
                compressedSignal(i) = intervalVector(j,1) + intervalVector(j,2);
                break;
            end
        end
    end
end