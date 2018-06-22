tmwaveform = csvread('OriginalSignal.csv');

HuffmanDynamicPolar(tmwaveform,32,40,1,false,false)

startA = 10;
startR = 6;
maxB = 40;
errorMax = 8;
error = zeros(maxB)+0;
avglen = zeros(maxB);
signalSize = zeros(maxB);
bitsMatrix = zeros(maxB);
exponent = zeros(maxB);
maxBits = ceil(log2(maxB^2))+1;
maxA = 1;
tmwaveform2 = normalization(tmwaveform);
huffman = false;
maxErrPlot = 10;

for i=startR:maxB
    for j=startA:maxB
        [error(i,j),avglen(i,j),signalSize(i,j)] = HuffmanDynamicPolar(tmwaveform,j,i,0,false,false);
        bitsMatrix(i,j) = i*j+1;
    end
end

intervalBits = 0:maxBits;
for i=1:maxBits
    logicalmatrix = (bitsMatrix > 2^intervalBits(i) & bitsMatrix <= 2^intervalBits(i+1));
    exponent = exponent + logicalmatrix .* 2^intervalBits(i+1);
end

dictUsage = bitsMatrix ./ exponent .*100;
wastedBits = exponent - bitsMatrix;
a = (error(startR:end,startA:end) <= maxErrPlot);
[rowI,columnI] = find(a == 1);
rowI = min(rowI)+startR;
columnI = min(columnI)+startA;
errorC = zeros(maxB) + max(max(error(rowI:end,columnI:end)));
errorC(rowI:end,columnI:end) = error(rowI:end,columnI:end);

imagesc(errorC);
colorbar
title('EVM(%) vs Num. Intervals Angles and Radius')
ylabel('Num. Angle')
xlabel('Num. Radius')
axis([columnI maxB rowI maxB])

if(huffman)
    avglenC = zeros(maxB) + max(max(avglen(rowI:end,columnI:end)));
    avglenC(rowI:end,columnI:end) = avglen(rowI:end,columnI:end);
    signalSizeC = zeros(maxB) + max(max(signalSize(rowI:end,columnI:end)));
    signalSizeC(rowI:end,columnI:end) = signalSize(rowI:end,columnI:end);
    
    figure
    imagesc(flip(avglenC));
    colorbar
    title('Average length vs Num. Intervals Angles and Radius')
    ylabel('Num. Angle')
    xlabel('Num. Radius')
    axis([columnI maxB 1 maxB-rowI])

    figure
    imagesc(flip(signalSizeC));
    colorbar
    title('Signal Size vs Num. Intervals Angles and Radius')
    ylabel('Num. Angle')
    xlabel('Num. Radius')
    axis([columnI maxB 1 maxB-rowI])
end

figure
contourf(dictUsage,1:10:100)
colorbar
title('Dictionary usage')
ylabel('Num. Angle')
xlabel('Num. Radius')
axis([startA maxB startR maxB])



minBits = min(bitsMatrix(errorC <= errorMax));
[row,column] = find(bitsMatrix == minBits & errorC <= errorMax);
[eee,aaa,sss] = HuffmanDynamicPolar(tmwaveform,column,row,0,false,false);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column),aaa,sss};