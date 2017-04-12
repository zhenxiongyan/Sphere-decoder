%%
clear all;
% Create an OFDM modulator and demodulator
ofdmMod = comm.OFDMModulator('FFTLength',128,'PilotInputPort',true,...
    'PilotCarrierIndices',cat(3,[7:6:120].',...
    [8:6:120].'),'InsertDCNull',true,...
    'NumTransmitAntennas',2);
ofdmDemod = comm.OFDMDemodulator(ofdmMod);
ofdmDemod.NumReceiveAntennas = 2;
% info(ofdmMod)
% disp(ofdmMod)
% info(ofdmDemod)
%%
% Show the resource mapping of pilot subcarriers for each transmit antenna.
%showResourceMapping(ofdmMod)
lteChannel = comm.LTEMIMOChannel(...
    'Profile',              'ETU 70HZ',...
    'AntennaConfiguration', '2x2',...
    'CorrelationLevel',     'Low',...
    'AntennaSelection',     'Off',...
    'RandomStream',         'mt19937ar with seed',...
    'Seed',                 1000,...
    'PathGainsOutputPort',  true,...
    'SampleRate', 1.92e6);
hAWGN = comm.AWGNChannel('NoiseMethod',...
    'Signal to noise ratio (SNR)','SNR',200);%%EbNo',20,'BitsPerSymbol',4)
%%
% Determine the dimensions of the OFDM modulator
ofdmModDim = info(ofdmMod);
numData = ofdmModDim.DataInputSize(1);   % Number of data subcarriers
numSym = ofdmModDim.DataInputSize(2);    % Number of OFDM symbols
numTxAnt = ofdmModDim.DataInputSize(3);  % Number of transmit antennas
%%
% Generate data symbols to fill 100 OFDM frames.
nframes = 100;modOrder = 2;
data = randi([0 1],nframes*numData*modOrder*2,1);
%% 16 QAM
hmod = comm.RectangularQAMModulator('ModulationOrder',2^modOrder,'BitInput',true,'NormalizationMethod','Average power');
Moddata = step(hmod,data(:));
%scatterplot(Moddata)
dataMod = reshape(Moddata,nframes*numData,numSym,numTxAnt);
%
%% OFDM Modulation
data_all = zeros(nframes*numData,numSym,numTxAnt);

for k = 1:nframes
    
    % Find row indices for kth OFDM frame
    indData = (k-1)*ofdmModDim.DataInputSize(1)+1:k*numData;
    
    % Generate random OFDM pilot symbols
    pilotData = complex(rand(ofdmModDim.PilotInputSize), ...
        rand(ofdmModDim.PilotInputSize));
    
    % Modulate  symbols using OFDM
    dataOFDM = step(ofdmMod,dataMod(indData,:,:),pilotData);
    
    %% LTE channel
    
    [LTEChanOut,LTEPathGain] = step(lteChannel,dataOFDM);
 
    %% AWGN Channel
    
    AWGNOut = step(hAWGN,LTEChanOut);
    
    %%
    % Demodulate OFDM data
    
    [receivedOFDMData,pilotOut] = step(ofdmDemod,AWGNOut);
    y=squeeze(receivedOFDMData);
    %%
    %%%   channel estimate
    dataCarrier_flag = true(ofdmMod.FFTLength,1);
    for tx = 1:2
        dataCarrier_flag(squeeze(ofdmMod.PilotCarrierIndices(:,1,tx))) = 0;
    end
    DC_indices = ofdmMod.FFTLength/2+1;
    GuardBand_indices = [1:ofdmMod.NumGuardBandCarriers(1) ofdmMod.FFTLength - ofdmMod.NumGuardBandCarriers(2)+1:ofdmMod.FFTLength].';
    dataCarrier_flag(DC_indices) = 0;
    dataCarrier_flag(GuardBand_indices) = 0;
    H_est = zeros(numData,2,2);
    for tx = 1:2
        for rx = 1:2
            % LS ????????
            H_LS = squeeze(pilotOut(:,1,tx,rx)./pilotData(:,1,tx));
            % ????
            p_idx = squeeze(ofdmMod.PilotCarrierIndices(:,1,tx));
            inter_p_idx = 1:ofdmMod.FFTLength;
            % ????
            H_LS_inter = interp1(p_idx,H_LS,inter_p_idx);
            % ????
            H_LS_inter(1:min(p_idx)-1) = H_LS_inter(min(p_idx));
            H_LS_inter(max(p_idx)+1:end) = H_LS_inter(max(p_idx));
            % IFFT
            h_ls = ifft(H_LS_inter);
            % ????
            [~,max_pos] = max(abs(h_ls));
            window = [mod(max_pos+20,128)+1 mod(max_pos-20,128)+1];
            h_ls(min(window):max(window)) = 0;
            % FFT
            H_LS = fft(h_ls);
            H_est(:,tx,rx) = H_LS(dataCarrier_flag);
        end
    end
    ChanG = H_est;
    
    %%
    
    %FSD detection
    numData = size(receivedOFDMData, 1);
    Variance = 1/hAWGN.SNR;
    r=zeros(numData,2);
    r_1=zeros(numData,2);
    X = complex(zeros(78,2));
    V = complex(zeros(78,2));
  
    for i = 1:numData
         h = squeeze(ChanG(i,:,:)).';        
        % FSD
        [Q,R] = qr(h);
        % the first stream
        X(i,2) = y(i,2)/R(2,2);
        % the second stream
        X(i,1) = (y(i,1) - R(1,2)*X(2,1))/R(1,1);
        if real(X(i,1))>0
            if imag(X(i,1))>0
                 X(i,1)=1+1j;
            else
                 X(i,1)=1-1j;
            end
        elseif real(X(i,1)) <0
            if imag(X(i,1))<0
                 X(i,1)=-1-1j;
            else
                 X(i,1)=-1+1j;
            end
        end   
         r(i,:)= X(i,:);
    end
    
    
    data_all(indData,:,:)= r;
    
end
%%



%% 16 demodulation
rxSig= data_all(:);
scatterplot(rxSig)
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',2^modOrder,'BitOutput',true);
dataout = step(hDemod,rxSig);
% err_bit = sum(data~=dataout)
% total_bit = length(dataout)
%
