function [result] = EVM(initialSignal,compressedSignal)

difference = abs(compressedSignal - initialSignal);




result = sqrt(sum(difference.^2)./...
    sum(abs(initialSignal).^2)) * 100;

end