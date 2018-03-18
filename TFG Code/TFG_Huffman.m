%% 
% Variables

roundFactor = 1;
%% 
% Plots w/o normalization
%%
plotSignal(tmwaveform);

% maxValueR = max(real(tmwaveform));
% maxValueI = min(real(tmwaveform));
%% 
% Normalization
%%
% maxValueNR = max([max(real(tmwaveform)),abs(min(real(tmwaveform)))]);
% maxValueNI = max([max(imag(tmwaveform)),abs(min(imag(tmwaveform)))]);
tmwaveform2 = normalization(tmwaveform);
%% 
% Plots regulated
%%
plotSignal(tmwaveform2);
%%
tmwaveform3 = splitIQ(tmwaveform2, roundFactor);
dictR = dictionaryEncoder(real(tmwaveform3),roundFactor);
dictI = dictionaryEncoder(imag(tmwaveform3),roundFactor);
plotSignal(tmwaveform3);
% %% 
% % Differencies
% %%
% different = zeros([length(tmwaveform3)-1 1]);
% for i = 1:length(tmwaveform3)-1
%     different(i) = real(tmwaveform3(i + 1)) + real(tmwaveform3 (i));
% end
% histogram(different);
% uniDiff = unique(different);
% %% 
% % Extra (rounded to 1)
% %%
% tmwaveform4 = round(tmwaveform2,1);
% numbins2 = min(real(tmwaveform4)): .1 : max(real(tmwaveform4));
% histogram(real(tmwaveform4),numbins2);
%% 
% Huffman coding
%%
%%[~,idx] = sort(Probalbility(:,2),'descend');
%%sortedProb = Probalbility(idx,:);
% infoRSended = huffmanenco(real(tmwaveform3'),dictR);
% infoISended = huffmanenco(imag(tmwaveform3'),dictI);
% % infoRReceived = huffmandeco(infoRSended',dictR);
% % infoIReceived = huffmandeco(infoISended',dictI);
% infoReceivedTotal = infoRSended + 1i * infoISended;
result = EVM(tmwaveform2,tmwaveform3);

