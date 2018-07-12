% tmwaveform = csvread('../../IFFT_OUTPUT.csv');
tmwaveform = csvread('OriginalSignal.csv');
% tmwaveform = csvread('../../DUC_OUTPUT.csv');
% tmwaveform = tmwaveform(1:30000);
% correctedR = real(tmwaveform) .* ~(real(tmwaveform) < -0.25) + -0.25 .* (real(tmwaveform) < -0.25);
% correctedI = imag(tmwaveform) .* ~(imag(tmwaveform) < -0.25) + -0.25 .* (imag(tmwaveform) < -0.25);
% corrected = correctedR + 1i *correctedI;
startVal = 10;
maxVal = 120;
errormaxB = 8;
error = zeros(1,maxVal)+500;
avglen = zeros(1,maxVal);
signalSize = zeros(1,maxVal);
bitsMatrix = zeros(1,maxVal);
huffman = true;
fit = false;

% stdSignalC = bitStd(corrected);
% stdSignalC = normalization(stdSignalC);

% stdSignal = bitStd(tmwaveform);
stdSignal = normalization(tmwaveform);
% plotSignal(stdSignal)

for i=startVal:maxVal
    [error(i),avglen(i),signalSize(i)] = HuffmanSplit(stdSignal,i,false,huffman);
    bitsMatrix(i) = i;
end
errorA = error;
avglenA = avglen;
signalSizeA = signalSize;

[dictUsage,wastedBits,exponent] = plot1DResults(error,bitsMatrix,avglen,signalSize,startVal,maxVal,huffman);

minBits = min(bitsMatrix(error <= errormaxB));
[row,column] = find(bitsMatrix == minBits & error <= errormaxB);
[eee,aaa,sss] = HuffmanSplit(stdSignal,minBits,true,true);
bestConf = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
    eee,log2(exponent(row,column)),bitsMatrix(row,column),...
    wastedBits(row,column),dictUsage(row,column),aaa,sss};

%% Fit
if(fit)
    [eee,aaa,sss] = HuffmanFitSplit(stdSignal,minBits,true,true);
    bestConf2 = {'Error','Num. Bits','Num Values','Wasted Values','DictUsage','Avg. len','Size Signal';...
        eee,log2(exponent(row,column)),bitsMatrix(row,column),...
        wastedBits(row,column),dictUsage(row,column),aaa,sss};
end