function [intervalVariable] = intervalVariable(interval,min,max)

intervalVariable = [min + interval/2:interval:max-interval/2]';
for i=1:length(intervalVariable)-1
    intervalVariable(i,2) = (intervalVariable(i+1,1) - intervalVariable(i,1))/2;
end