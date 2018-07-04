function [newData] = normalization(signalData)

maxValueNR = max([max(real(signalData)),abs(min(real(signalData)))]);
maxValueNI = max([max(imag(signalData)),abs(min(imag(signalData)))]);
if(maxValueNI == 0)
    newData = 1/maxValueNR * real(signalData);
elseif(maxValueNR == 0)
    newData = 0 + 1/maxValueNI * 1i * imag(signalData);
else
    newData = 1/maxValueNR * real(signalData) + 1/maxValueNI * 1i * imag(signalData);
end

end