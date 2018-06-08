% profile on
tmwaveform = csvread('OriginalSignal.csv');
startB = 350;
maxB = 400;
errormaxB = 8;
maxV = 20;
error = zeros(maxB,maxV)+500;
avglen = zeros(maxB,maxV);
signalSize = zeros(maxB,maxV);
bitsMatrix = zeros(maxB,maxV);
exponent = zeros(maxB,maxV);
maxBBits = ceil(log2(maxB));
huffman = false;
startVueltas = 5;


for i=startB:maxB
    for j=startVueltas:maxV
        [error(i,j),avglen(i,j),signalSize(i,j)] = HuffmanSpiral(tmwaveform,i,j,false,huffman);
        bitsMatrix(i,j) = i;
    end    
end

% profile viewer

intervalBits = 0:maxBBits;
for i=1:maxBBits
    logicalmatrix = (bitsMatrix > 2^intervalBits(i) & bitsMatrix <= 2^intervalBits(i+1));
    exponent = exponent + logicalmatrix .* 2^intervalBits(i+1);
end

dictUsage = bitsMatrix ./ exponent .*100;
wastedBits = exponent - bitsMatrix;
a = (error(startB:end,startVueltas:end) <= 15);
[rowI,columnI] = find(a == 1);
rowI = min(rowI)+startB;
columnI = min(columnI)+startVueltas;
errorC = zeros(maxB,maxV) + max(max(error(rowI:end,columnI:end)));
errorC(rowI:end,columnI:end) = error(rowI:end,columnI:end);

figure
imagesc(flip(errorC));
colorbar
title('EVM(%) vs Num. Values and Spin')
ylabel('Num. Values')
xlabel('Num. Spins')
axis([columnI maxV 1 maxB-rowI])

if(huffman)
    avglenC = zeros(maxB,maxV) + max(max(avglen(rowI:end,columnI:end)));
    avglenC(rowI:end,columnI:end) = avglen(rowI:end,columnI:end);
    signalSizeC = zeros(maxB,maxV) + max(max(signalSize(rowI:end,columnI:end)));
    signalSizeC(rowI:end,columnI:end) = signalSize(rowI:end,columnI:end);
    
    figure
    imagesc(flip(avglenC));
    colorbar
    title('Average length vs Num. Values and Spin')
    ylabel('Num. Values')
    xlabel('Num. Spins')
    axis([columnI maxV 1 maxB-rowI])

    figure
    imagesc(flip(signalSizeC));
    colorbar
    title('Signal Size vs Num. Values and Spin')
    ylabel('Num. Values')
    xlabel('Num. Spins')
    axis([columnI maxV 1 maxB-rowI])
end

figure
contourf(dictUsage,1:10:100)
colorbar
title('Dictionary usage')
ylabel('Num. Values')
xlabel('Num. Spins')
axis([startVueltas maxV startB maxB])

minBits = min(bitsMatrix(errorC <= errormaxB));
[row,column] = find(bitsMatrix == minBits & errorC <= errormaxB);
[eee,aaa,sss] = HuffmanSpiral(tmwaveform,row,column,true,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column),aaa,sss};

%% Quant fit
[eee,aaa,sss,newLen] = HuffmanFitSpiral(tmwaveform,row,column,true,true);
bestConf2 = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(newLen,column)),newLen,...
    wastedBits(newLen,column),dictUsage(newLen,column),aaa,sss};

% profile viewer