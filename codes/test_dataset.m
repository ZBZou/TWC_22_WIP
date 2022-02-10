%% 
%1:2048 data bits
%2049 M-QAM
%2050 number of Tx
%2051 Channel model
%2052 SNR
%2053 number of Rx
%2054:2565 received symbol
C = readtable('16QAM_TX4_RX4.csv'); 
%%
data_bit = C(1,1:2048);
data_bit = data_bit{:,:}';
rx_symb = C(1,2054:2565);
rx_symb = rx_symb{:,:}';

tx_symb = qammod(data_bit,16,'InputType','bit','UnitAveragePower',true);
%%
cd = comm.ConstellationDiagram('ShowReferenceConstellation',false);
cd(rx_symb)
