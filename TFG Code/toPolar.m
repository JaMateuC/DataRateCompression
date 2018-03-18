function polarSignal = toPolar(signal)

polarSignal = [abs(signal),atan2(imag(signal),real(signal))];
for i=1:length(polarSignal)
    if(polarSignal(i,2) < 0)
        polarSignal(i,2) = 2*pi + polarSignal(i,2);
    end
end