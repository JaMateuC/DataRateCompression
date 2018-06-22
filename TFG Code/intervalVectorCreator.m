function [intervalVector] = intervalVectorCreator(vector1,vector2)

intervalVector = ones(length(vector2)*length(vector1),2);
for i=1:length(vector2)
    
   intervalVector((i-1)*length(vector1)+1:(i)*length(vector1),1) = vector2(i)-intervalRadius/2;
    
end
for i=1:length(vector1)
    
   intervalVector(i:length(vector1):length(intervalVector),2) = vector1(i)-intervalAngle/2;
    
end