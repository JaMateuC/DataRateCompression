function [dictUsage,wastedBits,exponent] = plot1DResults(error,bitsMatrix,avglen,signalSize,startVal,maxVal,huffman)

exponent = zeros(1,maxVal);
maxBBits = ceil(log2(maxVal));

intervalBits = 0:maxBBits;
for i=1:maxBBits
    logicalmatrix = (bitsMatrix > 2^intervalBits(i) & bitsMatrix <= 2^intervalBits(i+1));
    exponent = exponent + logicalmatrix .* 2^intervalBits(i+1);
end

dictUsage = bitsMatrix ./ exponent .*100;
wastedBits = exponent - bitsMatrix;

figure
plot(error);
title('EVM(%) vs num. Values')
ylabel('EVM(%)')
xlabel('Num. Values')
axis([startVal maxVal 0 max(error(startVal:end))])

if(huffman)
    figure
    plot(avglen)
    title('Average length vs num. Values')
    xlabel('Intervals')
    ylabel('size (bits)')
    axis([startVal maxVal 0 max(avglen(startVal:end))])

    figure
    plot(signalSize)
    title('Signal size vs num. Values')
    xlabel('Intervals')
    ylabel('size(bits)')
    axis([startVal maxVal min(signalSize(startVal:end)) max(signalSize(startVal:end))])
end

figure
plot(dictUsage);
title('dictionary usage(%) vs num. Values')
ylabel('usage(%)')
xlabel('Num. Values')
axis([startVal maxVal min(dictUsage) 100])