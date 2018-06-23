tmwaveform = csvread('OriginalSignal.csv');
startVal = 350;
maxVal = 400;
errormaxB = 8;
maxSpin = 20;
error = zeros(maxVal,maxSpin)+500;
avglen = zeros(maxVal,maxSpin);
signalSize = zeros(maxVal,maxSpin);
bitsMatrix = zeros(maxVal,maxSpin);
huffman = true;
startSpins = 5;

tmwaveform2 = normalization(tmwaveform);
stdSignal = bitStd(tmwaveform2);

for i=startVal:maxVal
    for j=startSpins:maxSpin
        [error(i,j),avglen(i,j),signalSize(i,j)] = HuffmanSpiral(stdSignal,i,j,false,huffman);
        bitsMatrix(i,j) = i;
    end    
end

[errorC,dictUsage,wastedBits,exponent] = plotResults(error,bitsMatrix,avglen,signalSize,startVal,startSpins,maxVal,maxSpin,huffman,'spiral');

minBits = min(bitsMatrix(errorC <= errormaxB));
[row,column] = find(bitsMatrix == minBits & errorC <= errormaxB);
[eee,aaa,sss] = HuffmanSpiral(stdSignal,row,column,true,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column),aaa,sss};

%% Quant fit
[eee,aaa,sss,newLen] = HuffmanFitSpiral(stdSignal,row,column,true,true);
bestConf2 = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(newLen,column)),newLen,...
    wastedBits(newLen,column),dictUsage(newLen,column),aaa,sss};