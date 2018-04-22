% profile on
tmwaveform = csvread('OriginalSignal.csv');
maxB = 40;
errorMax = 8;
error = zeros(maxB);
avglen = zeros(maxB);
signalSize = zeros(maxB);
bitsMatrix = zeros(maxB);
exponent = zeros(maxB);
maxBits = ceil(log2(maxB^2));
startP = 10;
startQ = 10;
huffman = true;
maxErrPlot = 15;

for i=startP:maxB
    for j=startQ:maxB
        [error(i,j),avglen(i,j),signalSize(i,j)] = HuffmanIQTogether(tmwaveform,i,j,false,huffman);
        bitsMatrix(i,j) = i*j;
    end
end

intervalBits = 0:maxBits;
for i=1:maxBits
    logicalmatrix = (bitsMatrix > 2^intervalBits(i) & bitsMatrix <= 2^intervalBits(i+1));
    exponent = exponent + logicalmatrix .* 2^intervalBits(i+1);
end

dictUsage = bitsMatrix ./ exponent .*100;
wastedBits = exponent - bitsMatrix;
a = (error(startP:end,startQ:end) <= 15);
[~,rowI] = max(a);
[~,columnI] = max(a,[],2);
rowI = min(rowI)+startP-1;
columnI = min(columnI)+startQ-1;
errorC = zeros(maxB) + max(max(error(rowI:end,columnI:end)));
errorC(rowI:end,columnI:end) = error(rowI:end,columnI:end);

imagesc(flip(errorC));
colorbar
title('EVM(%) vs Num. Intervals Phase and Quadrature')
ylabel('Num. Phase')
xlabel('Num. Quadrature')
axis([columnI maxB 1 maxB-rowI])

if(huffman)
    avglenC = zeros(maxB) + max(max(avglen(rowI:end,columnI:end)));
    avglenC(rowI:end,columnI:end) = avglen(rowI:end,columnI:end);
    signalSizeC = zeros(maxB) + max(max(signalSize(rowI:end,columnI:end)));
    signalSizeC(rowI:end,columnI:end) = signalSize(rowI:end,columnI:end);
    
    figure
    imagesc(flip(avglenC));
    colorbar
    title('Average length vs Num. Intervals Phase and Quadrature')
    ylabel('Num. Phase')
    xlabel('Num. Quadrature')
    axis([columnI maxB 1 maxB-rowI])

    figure
    imagesc(flip(signalSizeC));
    colorbar
    title('Signal Size vs Num. Intervals Phase and Quadrature')
    ylabel('Num. Phase')
    xlabel('Num. Quadrature')
    axis([columnI maxB 1 maxB-rowI])
end

figure
contourf(dictUsage,1:10:100)
colorbar
title('Dictionary usage')
ylabel('Num. Phase')
xlabel('Num. Quadrature')
axis([startQ maxB startP maxB])

minBits = min(bitsMatrix(errorC <= errorMax));
[row,column] = find(bitsMatrix == minBits & errorC <= errorMax);
[eee,aaa,sss] = HuffmanIQTogether(tmwaveform,row,column,true,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column),aaa,sss};

% profile viewer