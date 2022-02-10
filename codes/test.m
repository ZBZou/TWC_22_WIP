%%WINNER II test

s = rng(10);                           % Set RNG seed for repeatability

% Transmission parameters
chanBW      = 'CBW80';               % Channel bandwidth
numUsers    = 3;                     % Number of users
numSTSVec   = [3 1 4];               % Number of streams per user
userPos     = [0 1 2];               % User positions
mcsVec      = [4 6 6];               % MCS per user: 16QAM, 64QAM, 256QAM
apepVec     = [355 180 600];         % Payload per user, in bytes
chCodingVec = {'BCC','LDPC','LDPC'}; % Channel coding per user

% Precoding and equalization parameters
precodingType = 'ZF';                % Precoding type; ZF or MMSE
snr           = 40;                  % SNR in dB
eqMethod      = 'ZF';                % Equalization method

% Create the multi-user VHT format configuration object
numTx = sum(numSTSVec);
cfgVHTMU = wlanVHTConfig('ChannelBandwidth',chanBW, ...
    'NumUsers',numUsers, ...
    'NumTransmitAntennas',numTx, ...
    'GroupID',2, ...
    'NumSpaceTimeStreams',numSTSVec,...
    'UserPositions',userPos, ...
    'MCS',mcsVec, ...
    'APEPLength',apepVec, ...
    'ChannelCoding',chCodingVec);

%%
% The number of transmit antennas is set to be the sum total of all the
% used space-time streams. This implies no space-time block coding (STBC)
% or spatial expansion is employed for the transmission.

%% Sounding (NDP) Configuration
%
% For precoding, channel sounding is first used to determine the channel
% experienced by the users (receivers). This channel state information is
% sent back to the transmitter, for it to be used for subsequent data
% transmission. It is assumed that the channel varies slowly over the two
% transmissions. For multi-user transmissions, the same NDP (Null Data
% Packet) is transmitted to each of the scheduled users [ <#18 2> ].

% VHT sounding (NDP) configuration, for same number of streams
cfgVHTNDP = wlanVHTConfig('ChannelBandwidth',chanBW, ...
    'NumUsers',1, ...
    'NumTransmitAntennas',numTx, ...
    'GroupID',0, ...
    'NumSpaceTimeStreams',sum(numSTSVec),...
    'MCS',0, ...
    'APEPLength',0);

%%
% The number of streams specified is the sum total of all space-time
% streams used. This allows the complete channel to be sounded.

% Generate the null data packet, with no data
txNDPSig = wlanWaveformGenerator([],cfgVHTNDP);
NPDSigLen = size(txNDPSig, 1);

%% WINNER II Channel for Indoor Office (A1) Scenario
%
% In this example, one |comm.WINNER2Channel| System object in the WINNER II
% Channel Model for Communications Toolbox(TM) is set up to simulate the
% three channels to different users. The indoor office (A1)
% non-line-of-sight (NLOS) scenario is configured for each user. With a
% fixed power delay profile, each user experiences a 16-path fading channel
% with the largest delay of 175 us. Each user is also assigned a low
% mobility as appropriate for 802.11ac.
%
% The access point employs a uniform circular array (UCA) with a radius of
% 20cm. Each user employs a uniform linear array (ULA) with 5cm spacing
% between elements. It is also assumed that each user's number of receive
% antennas is equal to the number of space-time streams allocated to them.

% Set up layout parameters for WINNER II channel
AA = winner2.AntennaArray('UCA',numTx,0.2);
for i = 1:numUsers
    AA(i+1) = winner2.AntennaArray('ULA',numSTSVec(i),0.05);
end
STAIdx   = 2:(numUsers+1);
APIdx   = {1};
rndSeed = 12;
cfgLayout = winner2.layoutparset(STAIdx,APIdx,numUsers,AA,[],rndSeed);
cfgLayout.Pairing = [ones(1,numUsers);2:(numUsers+1)]; % One access point to all users
cfgLayout.ScenarioVector = ones(1,numUsers);           % A1 scenario for all links
cfgLayout.PropagConditionVector = zeros(1,numUsers);  % NLOS
for i = 1:numUsers % Randomly set velocity for each user
    v = rand(3,1) - 0.5;
    cfgLayout.Stations(i+1).Velocity = v/norm(v,'fro');
end

% Set up model parameters for WINNER II channel
cfgModel = winner2.wimparset;
cfgModel.FixedPdpUsed       = 'yes';
cfgModel.FixedAnglesUsed    = 'yes';
cfgModel.IntraClusterDsUsed = 'no';
cfgModel.RandomSeed         = 111;    % Repeatability

% The maximum velocity for the 3 users is 1m/s. Set up the SampleDensity
% field to ensure that the sample rate matches the channel bandwidth.
maxMSVelocity = max(cell2mat(cellfun(@(x) norm(x,'fro'), ...
    {cfgLayout.Stations.Velocity},'UniformOutput',false)));
cfgModel.UniformTimeSampling = 'yes';
cfgModel.SampleDensity = round(physconst('LightSpeed')/ ...
    cfgModel.CenterFrequency/2/(maxMSVelocity/wlanSampleRate(cfgVHTMU)));

% Create the WINNER II channel System object
WINNERChan = comm.WINNER2Channel(cfgModel,cfgLayout);

% Call the info method to check some derived channel parameters
chanInfo = info(WINNERChan)

%%
% The channel filtering delay for each user is stored to account for its
% compensation at the receiver. In practice, symbol timing estimation would
% be used. At transmitter, an extra ten all-zero samples are appended to
% account for channel filter delay.

chanDelay   = chanInfo.ChannelFilterDelay;
numPadZeros = 10;

% Set ModelConfig.NumTimeSamples to match the length of the input signal to
% avoid warning
WINNERChan.ModelConfig.NumTimeSamples = NPDSigLen + numPadZeros;

% Sound the WINNER II channel for all users
chanOutNDP = WINNERChan([txNDPSig;zeros(numPadZeros,numTx)]);

% Add AWGN
rxNDPSig = cellfun(@awgn,chanOutNDP, ...
    num2cell(snr*ones(numUsers,1)),'UniformOutput',false);

%% Channel State Information Feedback
%
% Each user estimates its own channel using the received NDP signal and
% computes the channel state information that it can send back to the
% transmitter. This example uses the singular value decomposition of the
% channel seen by each user to compute the CSI feedback.

mat = cell(numUsers,1);
for uIdx = 1:numUsers
    % Compute the feedback matrix based on received signal per user
    mat{uIdx} = vhtCSIFeedback(rxNDPSig{uIdx}(chanDelay(uIdx)+1:end,:), ...
        cfgVHTNDP,uIdx,numSTSVec);
end

%%
% Assuming perfect feedback, with no compression or quantization loss of
% the CSI, the transmitter computes the steering matrix for the data
% transmission using either Zero-Forcing or Minimum-Mean-Square-Error
% (MMSE) based precoding techniques. Both methods attempt to cancel out the
% intra-stream interference for the user of interest and interference due
% to other users. The MMSE-based approach avoids the noise enhancement
% inherent in the zero-forcing technique. As a result, it performs better
% at low SNRs.

% Pack the per user CSI into a matrix
numST = length(mat{1});         % Number of subcarriers
steeringMatrix = zeros(numST,sum(numSTSVec),sum(numSTSVec));
%   Nst-by-Nt-by-Nsts
for uIdx = 1:numUsers
    stsIdx = sum(numSTSVec(1:uIdx-1))+(1:numSTSVec(uIdx));
    steeringMatrix(:,:,stsIdx) = mat{uIdx};     % Nst-by-Nt-by-Nsts
end

% Zero-forcing or MMSE precoding solution
if strcmp(precodingType, 'ZF')
    delta = 0; % Zero-forcing
else
    delta = (numTx/(10^(snr/10))) * eye(numTx); % MMSE
end
for i = 1:numST
    % Channel inversion precoding
    h = squeeze(steeringMatrix(i,:,:));
    steeringMatrix(i,:,:) = h/(h'*h + delta);
end

% Set the spatial mapping based on the steering matrix
cfgVHTMU.SpatialMapping = 'Custom';
cfgVHTMU.SpatialMappingMatrix = permute(steeringMatrix,[1 3 2]); %sub-by-Ns

%% Data Transmission
%
% Random bits are used as the payload for the individual users. A cell
% array is used to hold the data bits for each user, |txDataBits|. For a
% multi-user transmission the individual user payloads are padded such that
% the transmission duration is the same for all users. This padding process
% is described in Section 9.12.6 of [ <#18 1> ]. In this example for
% simplicity the payload is padded with zeros to create a PSDU for each
% user.

% Create data sequences, one for each user
txDataBits = cell(numUsers,1);
psduDataBits = cell(numUsers,1);
for uIdx = 1:numUsers
    % Generate payload for each user
    txDataBits{uIdx} = randi([0 1],cfgVHTMU.APEPLength(uIdx)*8,1,'int8');
    
    % Pad payload with zeros to form a PSDU
    psduDataBits{uIdx} = [txDataBits{uIdx}; ...
        zeros((cfgVHTMU.PSDULength(uIdx)-cfgVHTMU.APEPLength(uIdx))*8,1,'int8')];
end

%%
% Using the format configuration, |cfgVHTMU|, with the steering matrix, to
% generate the multi-user VHT waveform.

txSig = wlanWaveformGenerator(psduDataBits,cfgVHTMU);
%% %%%%%%%%%%%%%%%%%%%%%%%% HOGMT-Precoding %%%%%%%%%%%%%%%%%%%%

txDatasymbol = cell(numUsers,1);
txpadsymbol = cell(numUsers,1);
tx_s = zeros(sum(numSTSVec),numST);
for uIdx = 1:numUsers   
   txDatasymbol{uIdx} = qammod(txDataBits{uIdx},2^mcsVec(uIdx),...
      'InputType','bit','UnitAveragePower',true);
   Pad = numSTSVec(uIdx)*numST - length(txDatasymbol{uIdx});
   txpadsymbol{uIdx} = [txDatasymbol{uIdx};zeros(Pad,1)];
   if uIdx == 1
      tx_s(1:numSTSVec(uIdx),:) = reshape(txpadsymbol{uIdx},numSTSVec(uIdx),[]);
   else 
      tx_s(sum(numSTSVec(1:uIdx-1))+1:sum(numSTSVec(1:uIdx)),:) = reshape(txpadsymbol{uIdx},numSTSVec(uIdx),[]);
   end
end

h = zeros(sum(numSTSVec),sum(numSTSVec),numST);

for uIdx = 1:numUsers
    stsIdx = sum(numSTSVec(1:uIdx-1))+(1:numSTSVec(uIdx));
    h(stsIdx,:,:) = permute(mat{uIdx},[3,2,1]);   
end

X = zeros(sum(numSTSVec),numST);
ep = 0.0001;
for STid = 1:numST
   X(:,STid) = HOGMT(tx_s(:,STid), h(:,:,STid), ep);
end

%% Rx
for STid = 1:numST
   Rx(:,STid) = h(:,:,STid)*X(:,STid);
end

Rx_u = cell(numUsers,1);
for uIdx = 1:numUsers
   stsIdx = sum(numSTSVec(1:uIdx-1))+(1:numSTSVec(uIdx));
   temp = reshape(Rx(stsIdx,:),[],1);
   Rx_u{uIdx} = temp;
end

% Add AWGN
Rx_u = cellfun(@awgn,Rx_u, ...
    num2cell(snr*ones(numUsers,1)),'UniformOutput',false);

% demod
hat_s = cell(numUsers,1);

for uIdx = 1:numUsers   
   hat_s{uIdx} = qamdemod(Rx_u{uIdx},2^mcsVec(uIdx),...
      'OutputType','bit','UnitAveragePower',true);
end

%% Plt
for uIdx = 1:numUsers
    % User space-time streams
    stsU = numSTSVec(uIdx);
    s_u = reshape(Rx_u{uIdx}, numSTSVec(uIdx),[]);
    s_u(find(abs(s_u)<0.1)) = [];
    scaler(uIdx) = ceil(max(abs([real(s_u(:)); imag(s_u(:))])));
    for i = 1:stsU
        subplot(numUsers,max(numSTSVec),(uIdx-1)*max(numSTSVec)+i);
        plot(s_u,'.');
        axis square
        spAxes(sum([0 numSTSVec(1:(uIdx-1))])+i) = gca; % Store axes handle
        title(['User ' num2str(uIdx) ', Stream ' num2str(i)]);
        grid on;
    end
end
% Scale axes for all subplots and scale figure
for i = 1:numel(spAxes)
    xlim(spAxes(i),[-max(scaler) max(scaler)]);
    ylim(spAxes(i),[-max(scaler) max(scaler)]);
end

%% Ber
ber = inf(1, numUsers);
for uIdx = 1:numUsers
    idx = (1:cfgVHTMU.APEPLength(uIdx)*8).';
    [~,ber(uIdx)] = biterr(txDataBits{uIdx}(idx),hat_s{uIdx}(idx));
    disp(['HOGMT-Preocoding BER for User ' num2str(uIdx) ': ' num2str(ber(uIdx))]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The WINNER II channel object does not allow the input signal size to
% change once locked, so we have to call the release method before passing
% the waveform through it. In addition, as we restart the channel, we want
% it to re-process the NDP before the waveform so as to accurately mimic
% the channel continuity. Only the waveform portion of the channel's output
% is extracted for the subsequent processing of each user.

release(WINNERChan);

% Set ModelConfig.NumTimeSamples to match the total length of NDP plus
% waveform and padded zeros
WINNERChan.ModelConfig.NumTimeSamples = ...
    WINNERChan.ModelConfig.NumTimeSamples + length(txSig) + numPadZeros;

% Transmit through the WINNER II channel for all users, with 10 all-zero
% samples appended to account for channel filter delay
chanOut = WINNERChan([txNDPSig; zeros(numPadZeros,numTx); ...
    txSig; zeros(numPadZeros,numTx)]);

% Extract the waveform output for each user
chanOut = cellfun(@(x) x(NPDSigLen+numPadZeros+1:end,:),chanOut,'UniformOutput',false);

% Add AWGN
rxSig = cellfun(@awgn,chanOut, ...
    num2cell(snr*ones(numUsers,1)),'UniformOutput',false);

%% Data Recovery Per User
%
% The receive signals for each user are processed individually. The example
% assumes that there are no front-end impairments and that the transmit
% configuration is known by the receiver for simplicity.
%
% A user number specifies the user of interest being decoded for the
% transmission. This is also used to index into the vector properties of
% the configuration object that are user-specific.

% Get field indices from configuration, assumed known at receiver
ind = wlanFieldIndices(cfgVHTMU);

% Single-user receivers recover payload bits
rxDataBits = cell(numUsers,1);
% scaler = zeros(numUsers,1);
spAxes = gobjects(sum(numSTSVec),1);
hfig = figure('Name','Per-stream equalized symbol constellation');
for uIdx = 1:numUsers
    rxNSig = rxSig{uIdx}(chanDelay(uIdx)+1:end, :);
    
    % User space-time streams
    stsU = numSTSVec(uIdx);
    
    % Estimate noise power in VHT fields
    lltf = rxNSig(ind.LLTF(1):ind.LLTF(2),:);
    demodLLTF = wlanLLTFDemodulate(lltf,chanBW);
    nVar = helperNoiseEstimate(demodLLTF,chanBW,sum(numSTSVec));
    
    % Perform channel estimation based on VHT-LTF
    rxVHTLTF  = rxNSig(ind.VHTLTF(1):ind.VHTLTF(2),:);
    demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF,chanBW,numSTSVec);
    chanEst = wlanVHTLTFChannelEstimate(demodVHTLTF,chanBW,numSTSVec);
    
    % Recover information bits in VHT Data field
    rxVHTData = rxNSig(ind.VHTData(1):ind.VHTData(2),:);
    [rxDataBits{uIdx},~,eqsym] = wlanVHTDataRecover(rxVHTData, ...
        chanEst,nVar,cfgVHTMU,uIdx, ...
        'EqualizationMethod',eqMethod,'PilotPhaseTracking','None', ...
        'LDPCDecodingMethod','layered-bp','MaximumLDPCIterationCount',6);
    
    % Plot equalized symbols for all streams per user
%     scaler(uIdx) = ceil(max(abs([real(eqsym(:)); imag(eqsym(:))])));
    for i = 1:stsU
        subplot(numUsers,max(numSTSVec),(uIdx-1)*max(numSTSVec)+i);
        plot(reshape(eqsym(:,:,i),[],1),'.');
        axis square
        spAxes(sum([0 numSTSVec(1:(uIdx-1))])+i) = gca; % Store axes handle
        title(['User ' num2str(uIdx) ', Stream ' num2str(i)]);
        grid on;
    end
end

% Scale axes for all subplots and scale figure
for i = 1:numel(spAxes)
    xlim(spAxes(i),[-max(scaler) max(scaler)]);
    ylim(spAxes(i),[-max(scaler) max(scaler)]);
end
% pos = get(hfig,'Position');
% set(hfig,'Position',[pos(1)*0.7 pos(2)*0.7 1.3*pos(3) 1.3*pos(4)]);

%%
% Per-stream equalized symbol constellation plots validate the simulation
% parameters and convey the effectiveness of the technique. Note the
% discernible 16QAM, 64QAM and QPSK constellations per user as specified on
% the transmit end. Also observe the EVM degradation over the different
% streams for an individual user. This is a representative characteristic
% of the channel inversion technique.
%
% The recovered data bits are compared with the transmitted payload bits to
% determine the bit error rate.

% Compare recovered bits against per-user APEPLength information bits
ber = inf(1, numUsers);
for uIdx = 1:numUsers
    idx = (1:cfgVHTMU.APEPLength(uIdx)*8).';
    [~,ber(uIdx)] = biterr(txDataBits{uIdx}(idx),rxDataBits{uIdx}(idx));
    disp(['SOAT BER for User ' num2str(uIdx) ': ' num2str(ber(uIdx))]);
end
snr

rng(s); % Restore RNG state