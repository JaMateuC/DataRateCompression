function [intervalVector] = intervalVectorFun(vector1,vector2,interval1,interval2)

intervalVector = ones(length(vector1)*length(vector2),2);
for i=1:length(vector1)
    
   intervalVector((i-1)*length(vector2)+1:(i)*length(vector2),1) = vector1(i,1)-interval1/2;
    
end
for i=1:length(vector2)
    
   intervalVector(i:length(vector2):length(intervalVector),2) = vector2(i,1)-interval2/2;
    
end