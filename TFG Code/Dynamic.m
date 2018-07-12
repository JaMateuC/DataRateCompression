% tmwaveform = csvread('../../IFFT_OUTPUT.csv');
%  tmwaveform = csvread('OriginalSignal.csv');
tmwaveform = csvread('../../DUC_OUTPUT.csv');

startVal = 10;
maxVal = 120;
startTrueVal = 20;
errormaxB = 8;
error = zeros(1,maxVal)+500;
avglen = zeros(1,maxVal);
signalSize = zeros(1,maxVal);
bitsMatrix = zeros(1,maxVal);
huffman = true;
trueValueInterv = 0;

stdSignal = normalization(tmwaveform);
% stdSignal = bitStd(tmwaveform2);

for i=startVal:maxVal
    [error(i),avglen(i),signalSize(i)] = HuffmanDynamicSplit(stdSignal,i,trueValueInterv,false,huffman);
    bitsMatrix(:,i) = i;
end

[dictUsage,wastedBits,exponent] = plot1DResults(error,bitsMatrix,avglen,signalSize,startVal,maxVal,huffman);

minBits = min(bitsMatrix(error <= errormaxB));
[eee,aaa,sss] = HuffmanDynamicSplit(stdSignal,minBits,trueValueInterv,true,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(minBits)),bitsMatrix(minBits),...
    wastedBits(minBits),dictUsage(minBits),aaa,sss};