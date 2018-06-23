% Variables

tmwaveform = csvread('OriginalSignal.csv');
startAng = 10;
startRadius = 6;
maxVal = 40;
errorMax = 8;
error = zeros(maxVal);
signalSize = zeros(maxVal);
avglen = zeros(maxVal);
bitsMatrix = zeros(maxVal);
huffman = false;

% Algorithm

tmwaveform2 = normalization(tmwaveform);
stdSignal = bitStd(tmwaveform2);
polarSignal = toPolar(stdSignal);

for i=startRadius:maxVal
    for j=startAng:maxVal
        [error(i,j),avglen(i,j),signalSize(i,j)] = HuffmanPolar(polarSignal,tmwaveform2,i,j,false,huffman);
        bitsMatrix(i,j) = i*j+1;
    end
end

[errorC,dictUsage,wastedBits,exponent] = plotResults(error,bitsMatrix,avglen,signalSize,startRadius,startAng,maxVal,maxVal,huffman,'polar');

minBits = min(bitsMatrix(errorC <= errorMax));
[row,column] = find(bitsMatrix == minBits & errorC <= errorMax);
[eee,aaa,sss] = HuffmanPolar(polarSignal,tmwaveform2,row,column,true,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','Radius Values','Quadrature Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),row,column,dictUsage(row,column),aaa,sss};

%% Fit

[eee,aaa,sss] = HuffmanFitPolar(polarSignal,tmwaveform2,row,column,true,true);
bestConf2 = {'Error','Num. Bits','Num Values','Wasted Values','Radius Values','Quadrature Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column)-1,...
    wastedBits(row,column)+1,row,column,bitsMatrix(row,column)-1/exponent(row,column)*100,aaa,sss};