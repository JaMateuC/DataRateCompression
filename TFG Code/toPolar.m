function polarSignal = toPolar(signal)

polarSignal = [abs(signal),atan2(imag(signal),real(signal))];
polarSignal(:,2) = 2*pi*(polarSignal(:,2) < 0) + polarSignal(:,2);