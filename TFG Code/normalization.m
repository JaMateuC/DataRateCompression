function [newData] = normalization(signalData)

maxValueNR = max([max(real(signalData)),abs(min(real(signalData)))]);
maxValueNI = max([max(imag(signalData)),abs(min(imag(signalData)))]);
newData = 1/maxValueNR * real(signalData) + 1/maxValueNI * 1i * imag(signalData);

end