profile on
max = 40;
errorMax = 8;
error = zeros(max);
bitsMatrix = zeros(max);
exponent = zeros(max);
found = false;
maxBits = ceil(log2(max^2));
maxA = 1;
tmwaveform2 = normalization(tmwaveform);
polarSignal = toPolar(tmwaveform2);

for i=1:max
    for j=1:max
        error(i,j) = HuffmanPolar(polarSignal,tmwaveform2,i,j,false,maxA);
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

contourf(error,[0 8 9 10 11 12 13])
colorbar
title('Error vs Num. Intervals Angles and Radius')
ylabel('Num. Angles')
xlabel('Num. Radius')

figure
contourf(dictUsage,1:10:100)
colorbar
title('Dictionary usage')
ylabel('Num. Angles')
xlabel('Num. Radius')

minBits = min(bitsMatrix(error <= errorMax));
[row,column] = find(bitsMatrix == minBits & error <= errorMax);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage';...
    error(row,column),log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column)};

%HuffmanPolar(tmwaveform,row,column,true,maxA);

%% Next step

colMin = find(error(end,:)<=10, 1 );
rowMin = find(error(:,end)<=10, 1 );
errorN = error(rowMin:end,colMin:end);
exponentN = exponent(rowMin:end,colMin:end);
bitsMatrixN = bitsMatrix(rowMin:end,colMin:end);
wastedBitsN = wastedBits(rowMin:end,colMin:end);
dictUsageN = dictUsage(rowMin:end,colMin:end);
usefulMat = errorN > errorMax;

maxMin = [maxA,0.6];
intervals = [0.1,0.05,0.01,0.005,0.001];

for intv = 1:length(intervals)
    l =1;
    for maxA=maxMin(1):-intervals(intv):maxMin(2)
        for i=1:max-rowMin
            for j=1:max-colMin
                if(usefulMat(i,j))
                    errorN(i,j) = HuffmanPolar(polarSignal,tmwaveform2,i-1+rowMin,j-1+colMin,false,maxA);
                end
            end
        end
        minBits(l,1) = min(bitsMatrixN(errorN <= errorMax));
        minBits(l,2) = maxA;
        l=l+1;
    end
    [a,indx] = min(minBits(:,1));
    maxMin(1) = minBits(find(minBits(1:indx,1)>a,1,'last'),2);
    maxMin(2) = minBits(find(minBits(indx+1:end,1)>a,1)+indx,2);
    minBits = [];
end



% [row,column] = find(bitsMatrixN == minBits & errorN <= errorMax);
% bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage';...
%     errorN(row,column),log2(exponentN(row,column)),bitsMatrixN(row,column),...
%     wastedBitsN(row,column),dictUsageN(row,column)};
profile viewer