max = 40;
errorMax = 8;
error = zeros(max);
bitsMatrix = zeros(max);
exponent = zeros(max);
found = false;
maxBits = 0;

while(~found)
    
    maxBits = maxBits + 1;
    
    if(max^2 < 2^maxBits)
        found = true;
    end 
    
end

for i=1:max
    for j=1:max
        error(i,j) = HuffmanPolar(tmwaveform,i,j,false);
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
axis([20 max 7 max])
figure
contourf(dictUsage,1:10:100)

minBits = min(bitsMatrix(error <= errorMax));
[row,column] = find(bitsMatrix == minBits & error <= errorMax);
bestConf = {'Error','Exponent','NumBits','Wasted Bits','DictUsage';...
    error(row,column),exponent(row,column),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column)};

HuffmanPolar(tmwaveform,row,column,true);
