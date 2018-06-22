tmwaveform = csvread('OriginalSignal.csv');
tmwaveform2 = normalization(tmwaveform);

startVal = 10;
maxVal = 30;
startTrueVal = 20;
errormaxB = 8;
error = zeros(1,maxVal)+500;
avglen = zeros(1,maxVal);
signalSize = zeros(1,maxVal);
bitsMatrix = zeros(1,maxVal);
exponent = zeros(1,maxVal);
maxBBits = ceil(log2(maxVal));
huffman = true;
trueValue = 10;


for i=startVal:maxVal
    [error(i),avglen(i),signalSize(i)] = HuffmanDynamicSplit(tmwaveform,i,0,false,huffman);
    bitsMatrix(:,i) = i;
end

intervalBits = 0:maxBBits;
for i=1:maxBBits
    logicalmatrix = (bitsMatrix > 2^intervalBits(i) & bitsMatrix <= 2^intervalBits(i+1));
    exponent = exponent + logicalmatrix .* 2^intervalBits(i+1);
end

dictUsage = bitsMatrix ./ exponent .*100;
wastedBits = exponent - bitsMatrix;

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
    ylabel('bits')
    axis([startVal maxVal 0 max(avglen(startVal:end))])

    figure
    plot(signalSize)
    title('Signal size vs num. Values')
    xlabel('Intervals')
    ylabel('size(bits)')
    axis([startVal maxVal min(signalSize(startVal:end)) max(signalSize(startVal:end))])
end

minBits = min(bitsMatrix(error <= errormaxB));
[eee,aaa,sss] = HuffmanDynamicSplit(tmwaveform,minBits,0,true,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(minBits)),bitsMatrix(minBits),...
    wastedBits(minBits),dictUsage(minBits),aaa,sss};