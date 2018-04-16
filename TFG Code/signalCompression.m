function [compressedSignal] = signalCompression(signalInput,intervalVector,intervalVector2,max1V,min1V,max2V,min2V)

[minLen,indx] = min([size(intervalVector,1) size(intervalVector2,1)]);
switch indx
    case 1
        sVector = intervalVector;
        lVector = intervalVector2;
        signalInput1 = signalInput(:,1);
        signalInput2 = signalInput(:,2);
    case 2
        lVector = intervalVector;
        sVector = intervalVector2;
        signalInput1 = signalInput(:,2);
        signalInput2 = signalInput(:,1);
end
compressedSignalS = zeros(length(signalInput),1);
compressedSignalL = zeros(length(signalInput),1);
compressedSignal = zeros(length(signalInput),2);

for j=1:minLen-1
    compressedSignalS = (sVector(j,1) + sVector(j,2)).*(signalInput1 > sVector(j,1) & signalInput1 <= sVector(j+1,1)) +...
        ~(signalInput1 > sVector(j,1) & signalInput1 <= sVector(j+1,1)).*compressedSignalS;
    compressedSignalL = (lVector(j,1) + lVector(j,2)).*(signalInput2 > lVector(j,1) & signalInput2 <= lVector(j+1,1)) +...
        ~(signalInput2 > lVector(j,1) & signalInput2 <= lVector(j+1,1)).*compressedSignalL;
end
for j=minLen:size(lVector,1)-1
    compressedSignalL = (lVector(j,1) + lVector(j,2)).*(signalInput2 > lVector(j,1) & signalInput2 <= lVector(j+1,1)) +...
        ~(signalInput2 > lVector(j,1) & signalInput2 <= lVector(j+1,1)).*compressedSignalL;
end

switch indx
    case 1
        compressedSignal(:,1) = compressedSignalS;
        compressedSignal(:,2) = compressedSignalL;
    case 2
        compressedSignal(:,2) = compressedSignalS;
        compressedSignal(:,1) = compressedSignalL;
end

compressedSignal(:,1) = max1V.*(signalInput(:,1) > intervalVector(end,1)) + ~(signalInput(:,1) > intervalVector(end,1)).*compressedSignal(:,1);
compressedSignal(:,1) = min1V.*(signalInput(:,1) <= intervalVector(1,1)) + ~(signalInput(:,1) <= intervalVector(1,1)).*compressedSignal(:,1);
compressedSignal(:,2) = max2V.*(signalInput(:,2) > intervalVector2(end,1)) + ~(signalInput(:,2) > intervalVector2(end,1)).*compressedSignal(:,2);
compressedSignal(:,2) = min2V.*(signalInput(:,2) <= intervalVector2(1,1)) + ~(signalInput(:,2) <= intervalVector2(1,1)).*compressedSignal(:,2);