function [errorC,dictUsage,wastedBits,exponent] = plotResults(error,bitsMatrix,avglen,signalSize,startRow,startCol,maxValRow,maxValCol,huffman,type) 

exponent = zeros(maxValRow,maxValCol);
maxBits = ceil(log2(maxValRow^2))+1;
maxErrPlot = 15;

intervalBits = 0:maxBits;
for i=1:maxBits
    logicalmatrix = (bitsMatrix > 2^intervalBits(i) & bitsMatrix <= 2^intervalBits(i+1));
    exponent = exponent + logicalmatrix .* 2^intervalBits(i+1);
end

dictUsage = bitsMatrix ./ exponent .*100;
wastedBits = exponent - bitsMatrix;
a = (error(startRow:end,startCol:end) <= maxErrPlot);
[rowI,columnI] = find(a == 1);
rowI = min(rowI)+startRow;
columnI = min(columnI)+startCol;
errorC = zeros(maxValRow,maxValCol) + max(max(error(rowI:end,columnI:end)));
errorC(rowI:end,columnI:end) = error(rowI:end,columnI:end);

imagesc(errorC);
colorbar
switch type
    case 'polar'
        title('EVM(%) vs Num. Intervals Angles and Radius')
        ylabel('Num. Angle')
        xlabel('Num. Radius')
    case 'cartesian'
        title('EVM(%) vs Num. Intervals Phase and Quadrature')
        ylabel('Num. Phase')
        xlabel('Num. Quadrature')
    case 'spiral'
        title('EVM(%) vs Num. Values and Spin')
        ylabel('Num. Values')
        xlabel('Num. Spins')
end
axis([columnI maxValCol rowI maxValRow])

figure
switch type
    case 'polar'
        contourf(dictUsage,1:10:100)
        colorbar
        title('Dictionary usage')
        ylabel('Num. Angle')
        xlabel('Num. Radius')
        axis([startCol maxValCol startRow maxValRow])
    case 'cartesian'
        contourf(dictUsage,1:10:100)
        colorbar
        title('Dictionary usage')
        ylabel('Num. Phase')
        xlabel('Num. Quadrature')
        axis([startCol maxValCol startRow maxValRow])
    case 'spiral'
        plot(dictUsage(1:maxValRow,maxValCol))
        title('Dictionary usage')
        ylabel('Percentage')
        xlabel('Num. Values')
        axis([startRow maxValRow -inf inf])
end


if(huffman)
    avglenC = zeros(maxValRow,maxValCol) + max(max(avglen(rowI:end,columnI:end)));
    avglenC(rowI:end,columnI:end) = avglen(rowI:end,columnI:end);
    signalSizeC = zeros(maxValRow,maxValCol) + max(max(signalSize(rowI:end,columnI:end)));
    signalSizeC(rowI:end,columnI:end) = signalSize(rowI:end,columnI:end);
    
    figure
    switch type
    case 'polar'
        imagesc(avglenC);
        colorbar
        title('Average length vs Num. Angle and Radius')
        ylabel('Num. Angle')
        xlabel('Num. Radius')
        axis([columnI maxValCol rowI maxValRow])
    case 'cartesian'
        imagesc(avglenC);
        colorbar
        title('Average length vs Num. Phase and Quadrature')
        ylabel('Num. Phase')
        xlabel('Num. Quadrature')
        axis([columnI maxValCol rowI maxValRow])
    case 'spiral'
        plot(avglenC(1:maxValRow,maxValCol))
        title('Average length vs Num. Values')
        ylabel('Average length')
        xlabel('Num. Values')
        axis([startRow+1 maxValRow -inf inf])
    end

    figure
    switch type
    case 'polar'
        imagesc(signalSizeC);
        colorbar
        title('Signal Size (bits) vs Num. Angle and Radius')
        ylabel('Num. Angle')
        xlabel('Num. Radius')
        axis([columnI maxValCol rowI maxValRow])
    case 'cartesian'
        imagesc(signalSizeC);
        colorbar
        title('Signal Size (bits) vs Num. Phase and Quadrature')
        ylabel('Num. Phase')
        xlabel('Num. Quadrature')
        axis([columnI maxValCol rowI maxValRow])
    case 'spiral'
        plot(signalSizeC(1:maxValRow,maxValCol))
        title('Signal Size (bits) vs Num. Values')
        xlabel('Num. Values')
        ylabel('Signal Size')
        axis([startRow+1 maxValRow -inf inf])
    end
end