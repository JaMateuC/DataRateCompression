tmwaveform = csvread('OriginalSignal.csv');

maxVal = 40;
errorMax = 8;
error = zeros(maxVal);
avglen = zeros(maxVal);
signalSize = zeros(maxVal);
bitsMatrix = zeros(maxVal);
startP = 10;
startQ = 10;
huffman = false;

tmwaveform2 = normalization(tmwaveform);
stdSignal = bitStd(tmwaveform2);

for i=startP:maxVal
    for j=startQ:maxVal
        [error(i,j),avglen(i,j),signalSize(i,j)] = HuffmanIQTogether(stdSignal,i,j,false,huffman);
        bitsMatrix(i,j) = i*j;
    end
end

[errorC,dictUsage,wastedBits,exponent] = plotResults(error,bitsMatrix,avglen,signalSize,startP,startQ,maxVal,maxVal,huffman,'cartesian');

minBits = min(bitsMatrix(errorC <= errorMax));
[row,column] = find(bitsMatrix == minBits & errorC <= errorMax);
[eee,aaa,sss] = HuffmanIQTogether(stdSignal,row,column,false,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','Phase Values','Quadrature Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),row,column,dictUsage(row,column),aaa,sss};

%% Fit

[eee,aaa,sss] = HuffmanFitIQTogether(stdSignal,row,column,true,true);
bestConf2 = {'Error','Num. Bits','Num Values','Wasted Values','Phase Values','Quadrature Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),row,column,dictUsage(row,column),aaa,sss};
