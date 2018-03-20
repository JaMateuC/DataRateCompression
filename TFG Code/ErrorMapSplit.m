max = 512;
errorMax = 8;
error = zeros(1,max);
bitsMatrix = zeros(1,max);
exponent = zeros(1,max);
found = false;
maxBits = ceil(log2(max));


for i=1:max
    error(i) = HuffmanSplit(tmwaveform,i,false);
    bitsMatrix(i) = i;
end

intervalBits = 0:maxBits;
for i=1:maxBits
    logicalmatrix = (bitsMatrix > 2^intervalBits(i) & bitsMatrix <= 2^intervalBits(i+1));
    exponent = exponent + logicalmatrix .* 2^intervalBits(i+1);
end

dictUsage = bitsMatrix ./ exponent .*100;
wastedBits = exponent - bitsMatrix;

figure
plot(error)
title('EVM vs num. min. intervals')
xlabel('Intervals')
ylabel('EVM')
figure
plot(dictUsage)
title('Dictionary usage')
xlabel('Intervals')
ylabel('Percentage')

minBits = min(bitsMatrix(error <= errorMax));
[row,column] = find(bitsMatrix == minBits & error <= errorMax);
bestConf = {'Error','Exponent','NumBits','Wasted Bits','DictUsage';...
    error(row,column),exponent(row,column),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column)};

HuffmanSplit(tmwaveform,minBits,true);