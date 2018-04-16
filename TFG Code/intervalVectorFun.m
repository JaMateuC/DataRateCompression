function [intervalVector] = intervalVectorFun(vector1,vector2)

vector1(length(vector1) + 1,1) =  vector1(end,1) + vector1(end,2)*2;
vector1(end,2) =  vector1(length(vector1)-1,2);
vector2(length(vector2) + 1,1) =  vector2(end,1) + vector2(end,2)*2;
vector2(end,2) =  vector2(length(vector2)-1,2);

intervalVector = ones(length(vector1)*length(vector2),2);
for i=2:length(vector1)
    
   intervalVector((i-1)*length(vector2)+1:(i)*length(vector2),1) = vector1(i,1)-vector1(i,2);
    
end
for i=1:length(vector2)
    
   intervalVector(i:length(vector2):length(intervalVector),2) = vector2(i,1)-vector2(i,2);
    
end