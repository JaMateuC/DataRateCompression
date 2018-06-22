function [newData] = splitIQ(signalData, roundFactor)

newData = round(signalData,roundFactor);
interval = 1/10^roundFactor;
x = [0, length(newData)];

figure
subplot(2,1,1)
hold on
plot(real(newData));
for i = -1:interval: 1
    y = [i i];
    plot(x,y, 'Color','green')
end
title('In-Phase Signal')
xlim(x)
hold off

subplot(2,1,2)
hold on
plot(imag(newData));
title('In-Quadrature Signal')
xlim(x)
hold off

end