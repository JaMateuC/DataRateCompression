tmwaveform = csvread('OriginalSignal.csv');

startA = 10;
startR = 6;
maxVal = 40;
errorMax = 8;
error = zeros(maxVal);
avglen = zeros(maxVal);
signalSize = zeros(maxVal);
bitsMatrix = zeros(maxVal);
huffman = false;

tmwaveform2 = normalization(tmwaveform);
stdSignal = bitStd(tmwaveform2);

for i=startR:maxVal
    for j=startA:maxVal
        [error(i,j),avglen(i,j),signalSize(i,j)] = HuffmanDynamicPolar(stdSignal,j,i,0,false,false);
        bitsMatrix(i,j) = i*j+1;
    end
end

[errorC,dictUsage,wastedBits,exponent] = plotResults(error,bitsMatrix,avglen,signalSize,startR,startA,maxVal,maxVal,huffman,'polar');

minBits = min(bitsMatrix(errorC <= errorMax));
[row,column] = find(bitsMatrix == minBits & errorC <= errorMax);
[eee,aaa,sss] = HuffmanDynamicPolar(stdSignal,column,row,0,false,false);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column),aaa,sss};