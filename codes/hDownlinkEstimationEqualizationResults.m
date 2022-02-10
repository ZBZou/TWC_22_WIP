%hDownlinkEstimationEqualizationResults Plot the resource grids
%   hDownlinkEstimationEqualizationResults(RXGRID, EQGRID) plots the
%   received and equalized resource grids. This requires the receive grid
%   RXGRID and equalized grid EQGRID.

%   Copyright 2009-2014 The MathWorks, Inc.

function hDownlinkEstimationEqualizationResults(rxGrid, eqGrid)
    
    % Plot received grid error on logarithmic scale
    figure;
    dims = size(rxGrid);
    surf(20*log10(abs(rxGrid)));
    title('Received resource grid');
    ylabel('Subcarrier');
    xlabel('Symbol');
    zlabel('absolute value (dB)');
    axis([1 dims(2) 1 dims(1) -40 10]);

    % Plot equalized grid error on logarithmic scale
    figure
    surf(20*log10(abs(eqGrid)));
    title('Equalized resource grid');
    ylabel('Subcarrier');
    xlabel('Symbol');
    zlabel('absolute value (dB)');
    axis([1 dims(2) 1 dims(1) -40 10]);

end