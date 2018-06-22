function [stdSignal] = bitStd(signal)

numvalues = 2^16-1;

intv = [-1:2/numvalues:1];

realPart = real(signal);
imagPart = imag(signal);

realMin = abs(realPart(1:100) - intv);
imagMin = abs(imagPart(1:100) - intv);
[~,minIndR] = min(realMin,[],2);
[~,minIndI] = min(imagMin,[],2);
tmwaveformR = intv(minIndR)';
tmwaveformI = intv(minIndI)';

for i = 101:1:length(signal)
realMin = abs(realPart(i) - intv);
imagMin = abs(imagPart(i) - intv);
[~,minIndR] = min(realMin,[],2);
[~,minIndI] = min(imagMin,[],2);
tmwaveformR = [tmwaveformR;intv(minIndR)'];
tmwaveformI = [tmwaveformI;intv(minIndI)'];
end

stdSignal = tmwaveformR + 1i * tmwaveformI;

