function [interval] = intervalVariable(interval)

for i=1:length(interval)-1
    interval(i,2) = (interval(i+1,1) - interval(i,1))/2;
end