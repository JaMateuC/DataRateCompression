% profile on
tmwaveform = csvread('OriginalSignal.csv');
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
polarSignal = toPolar(tmwaveform2);
huffman = false;
maxErrPlot = 15;

for i=startR:maxB
    for j=startA:maxB
        [error(i,j),avglen(i,j),signalSize(i,j)] = HuffmanPolar(polarSignal,tmwaveform2,i,j,false,maxA,huffman);
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
a = (error(startR:end,startA:end) <= 15);
[rowI,columnI] = find(a == 1);
rowI = min(rowI)+startR;
columnI = min(columnI)+startA;
errorC = zeros(maxB) + max(max(error(rowI:end,columnI:end)));
errorC(rowI:end,columnI:end) = error(rowI:end,columnI:end);

imagesc(flip(errorC));
colorbar
title('EVM(%) vs Num. Intervals Angles and Radius')
ylabel('Num. Angle')
xlabel('Num. Radius')
axis([columnI maxB 1 maxB-rowI])

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
[eee,aaa,sss] = HuffmanPolar(polarSignal,tmwaveform2,row,column,true,maxA,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column),aaa,sss};

% profile viewer
%% Next step

% colMin = find(errorC(end,:)<=10, 1 );
% rowMin = find(errorC(:,end)<=10, 1 );
% errorN = errorC(rowMin:end,colMin:end);
% exponentN = exponent(rowMin:end,colMin:end);
% bitsMatrixN = bitsMatrix(rowMin:end,colMin:end);
% wastedBitsN = wastedBits(rowMin:end,colMin:end);
% dictUsageN = dictUsage(rowMin:end,colMin:end);
% avglenN = avglen(rowMin:end,colMin:end);
% signalSizeN = signalSize(rowMin:end,colMin:end);
% usefulMat = errorN > errorMax;
% 
% maxMin = [maxA,0.6];
% intervals = [0.1,0.05,0.01,0.005,0.001];
% 
% for intv = 1:length(intervals)
%     l =1;
%     minBits = zeros((maxMin(2)-maxMin(1))/intervals(intv),3);
%     for maxA=maxMin(1):-intervals(intv):maxMin(2)
%         for i=1:maxB-rowMin
%             for j=1:maxB-colMin
%                 if(usefulMat(i,j))
%                     [errorN(i,j),avglenN(i,j),signalSizeN(i,j)] = HuffmanPolar(polarSignal,tmwaveform2,i-1+rowMin,j-1+colMin,false,maxA,false);
%                 end
%             end
%         end
%         minBits(l,1) = min(bitsMatrixN(errorN <= errorMax));
%         minBits(l,2) = maxA;
%         minBits(l,3) = min(errorN(bitsMatrixN == minBits(l,1) & errorN <= errorMax));
%         [row,column] = find(bitsMatrixN == minBits(l,1) & errorN <= errorMax & errorN == minBits(l,3));
%         minBits(l,4) = row + rowMin;
%         minBits(l,5) = column + colMin;
%         l=l+1;
%     end
%     [a,indx] = min(minBits(:,1));
%     maxMin(1) = minBits(find(minBits(1:indx,1)>=a,1,'last'),2);
%     maxMin(2) = minBits(find(minBits(indx+1:end,1)>=a,1)+indx,2);
% end
% 
% mtxInt = minBits(minBits(:,1) == min(minBits(:,1)) ,:);
% [minE,indxInt] = min(mtxInt(:,3));
% minInt = mtxInt(indxInt,:);
% 
% [errorN,avglen,signalSize] = HuffmanPolar(polarSignal,tmwaveform2,minInt(4),minInt(5),false,minInt(2),true);
% 
% bestConf2 = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage',...
%     'Avg. len','Size Signal';errorN,log2(exponentN(row,column)),...
%     bitsMatrixN(row,column),wastedBitsN(row,column),...
%     dictUsageN(row,column),avglen,signalSize};
% 
% profile viewer