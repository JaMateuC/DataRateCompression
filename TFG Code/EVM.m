function [result] = EVM(initialSignal,compressedSignal,plots)

difference = abs(compressedSignal - initialSignal);

result = sqrt(sum(difference.^2)./...
    sum(abs(initialSignal).^2)) * 100;

if(plots)
    figure
    plot(difference)
    title('Error of each sample')
    xlabel('samples')
    ylabel('error')
end

end