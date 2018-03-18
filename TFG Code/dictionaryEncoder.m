function [dict] = dictionaryEncoder(signalIorQ,roundFactor)

rF = 1/(10^roundFactor);

numbins = -max(signalIorQ): rF : max(signalIorQ);
edges = -rF/2-max(signalIorQ): rF : max(signalIorQ)+rF/2;
figure
histogram(signalIorQ,numbins);
numCount = histcounts(signalIorQ,edges);
totalCount = sum(numCount);
Probalbility = [numbins',numCount'/totalCount];
dict = huffmandict(Probalbility(:,1),Probalbility(:,2));
for i=1:size(dict,1)
    dict{i,1} = round(dict{i,1},roundFactor);
end