file = 'IFFT_OUTPUT ';

fileID = fopen(['C:\Users\usuario\Desktop\'  file  '.txt'],'r');

formatSpec = '%f';
A = fscanf(fileID,formatSpec);

fclose(fileID);

signal = zeros(length(A)/2,2);

signal(:,1) = A(1:2:end);
signal(:,2) = A(2:2:end);

signalF = signal(:,1) + 1i* signal(:,2);

csvwrite(['C:\Users\usuario\Desktop\'  file  '.csv'],signalF)