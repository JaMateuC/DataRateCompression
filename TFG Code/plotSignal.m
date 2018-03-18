function [] = plotSignal(signalData)

figure
subplot(2,2,2)
plot(signalData,'x')
subplot(2,2,1)
histogram(imag(signalData));
view(-90, 90);
xlim([min(imag(signalData)) max(imag(signalData))]);
subplot(2,2,4)
histogram(real(signalData));
set(gca,'ydir','r');
xlim([min(imag(signalData)) max(imag(signalData))]);

end