close all
clear all
clc

Mtotal=4267;
Mdata=240;
Mdc=0;
Ncp=0;
Nsym = 600; 
fmezinos = 120e3;
fs = 512e6;
delka = 0.0200015625/4;
angles = 63;
smer = 1:angles;

%filename_tx = sprintf('data/tx_data/sig_tx_BFK_RX_state.mat');
filename_tx = sprintf('data/tx_data/sig_tx_final_09052025.mat');
load(filename_tx);

[M5_TX] = OFDM_demodulator(sig_tx, Mtotal, Mdata, Mdc, Ncp, Nsym);

% figure(1)
% imagesc((abs(M5_TX)));
% axis xy;
% xlabel('OFDM sym');
% ylabel('Subcar');
% title('SS burst, demodulated from RX sig_rx');
% colorbar;

K = 0;
a_TX = (floor(size(M5_TX, 1)/2+1)-240/2)-Mdc+K;
b_TX = (floor(size(M5_TX, 1)/2)+240/2)+Mdc+K;
M6_TX = fftshift(M5_TX, 1);

% figure(2)
% imagesc((abs(M6_TX)));
% axis xy;
% xlabel('OFDM sym');
% ylabel('Subcar');
% title('SS burst, demodulated from RX sig_rx and rearanged with fftshift');
% colormap jet;
% colorbar;

M7_TX = M6_TX(a_TX:b_TX,1:600);
M7_TX = [M7_TX(1:120,:); M7_TX(120+2*Mdc+1:end,:)];

% figure(3)
% imagesc((abs(M7_TX)));
% axis xy;
% xlabel('OFDM sym');
% ylabel('Subcar');
% colormap jet;
% colorbar;

M9_TX = M7_TX(:,[4:11, 16:23, 32:39, 44:51]);
M9_TX = M9_TX/max(max(abs(M9_TX)));
M9_TX = round(M9_TX,3);

figure(4)
imagesc((abs(M9_TX)));
axis xy;
set(gca, 'FontSize', 30);
xlabel('OFDM symbol','FontSize',35);
ylabel('Subcarrier frequency','FontSize',35);
title('8 Synchronization Signal Blocks','FontSize',35);
colormap jet;
clim([0 1])
colorbar;

for i = 31%1:angles
    disp("Nacitam soubor data_rx_" + i);
    if i > 1
        clear data sig_rx
    end
    filename_rx = sprintf('data/measure_15_4_2025/LOS/data_rx_%d.mat', i);
    %filename_rx = sprintf('data/measure_15_4_2025/NLOS/data_rx_%d.mat', i);
    load(filename_rx);
    sig_rx = data;

    %% time
    t = 0:1/fs:delka*4-(1/fs); % s
    t = t*1000; % ms

    %% display of received signal
    figure(5)
    %subplot(7,10,i)
    plot(t, abs(sig_rx))
    set(gca, 'FontSize', 30);
    xlabel('t [ms]', 'FontSize', 35);
    ylabel('Amplitude [-]', 'FontSize', 35);
    title('SS burst in time domain RX', 'FontSize', 35);
    xlim([0 20])

    %% compensation of f offset
    best_distances = inf(240,32);
    krok_fo = 1;
    for fo = -37e3:50:-33e3
        % frequency shift by x Hz
        % Parametres
        fshift = fo;      % frequency offset (Hz)
        fs = 512e6;        % sampling frequency
        T = 1/fs;          % period
        L = length(sig_rx);% number of samples

        % Block size
        block_size = 1e6;

        % Block processing
        for i_block = 1:block_size:L
            i_end = min(i_block + block_size - 1, L);
            n = (i_block-1):(i_end-1);                   % indices within the entire signal
            t_block = n * T;                             % time stamps
            sig_rx(i_block:i_end) = sig_rx(i_block:i_end) .* exp(1i * 2 * pi * fshift * t_block);
        end

        %% synchronization
        [korelace, lags] = xcorr(abs(sig_rx), abs(sig_tx));

        % figure(6);
        % %subplot(7,10,i)
        % plot(lags, korelace);
        % xlabel('lags');
        % ylabel('corelation');
        % title('Corelation of transmitted and recieved signal');

        %extract the original length for 5 ms from sig_rx
        [M, I] = max(korelace);
        k = 0;

        if lags(I)+k+length(sig_tx) > 4*length(sig_tx)
            I = I - length(sig_tx);
        elseif lags(I) < 0
            I = I + length(sig_tx);
        end

        sig_rx = [sig_rx(lags(I)+1+k:end)]; % cut out the synchronized window
        sig_rx = sig_rx/max(abs(sig_rx)); % normalization

        % figure(7)
        % plot(abs(sig_rx))
        % xlabel('t [ms]');
        % ylabel('amp');
        % title('Synchronized SS burst in time domain RX');

        %% synchronization comparison
        % figure(8);
        % %subplot(12,10,i)
        % plot(abs(sig_rx));
        % hold on;
        % plot(abs(sig_tx));
        % legend('RX', 'TX');
        % title('Time'); hold off;

        % figure(9);
        % %subplot(12,10,i)
        % plot(angle(sig_rx));
        % hold on;
        % plot(angle(sig_tx));
        % legend('RX phase', 'TX phase');
        % title('Phase'); hold off;

        %% demodulator
        [M5_RX] = OFDM_demodulator(sig_rx, Mtotal, Mdata, Mdc, Ncp, Nsym);

        K = 2; %constant for rearrangement in the spectrum
        a_RX = (floor(size(M5_RX, 1)/2+1)-240/2)-Mdc+K;
        b_RX = (floor(size(M5_RX, 1)/2)+240/2)+Mdc+K;

        %% fftshift - rotate the spectrum back to make it symmetrical around zero
        M6_RX = fftshift(M5_RX, 1);

        % figure(11)
        % %subplot(7,10,i)
        % imagesc((abs(M6_RX)));
        % axis xy;
        % xlabel('OFDM sym');
        % ylabel('Subcar');
        % title('SS burst, demodulated from RX sig_rx and rearanged with fftshift');
        % colormap jet;
        % colorbar;

        %% direction decoding
        M7_RX = M6_RX(a_RX:b_RX,1:600); %detail of the area where SSBurst is
        M7_RX = [M7_RX(1:120,:); M7_RX(120+2*Mdc+1:end,:)];

        % figure(12)
        % %subplot(7,10,i)
        % imagesc((abs(M7_RX)));
        % axis xy;
        % xlabel('OFDM sym');
        % ylabel('Subcar');
        % colormap jet;
        % colorbar;

        %% Reference Signal Received Power (RSRP)
        % stack 8 SSBlock symbols behind and calculate the RSRP for them
        M9_RX = M7_RX(:,[4:11, 16:23, 32:39, 44:51]);% M7_RX(:,5);

        % figure(13)
        % %subplot(7,10,i)
        % imagesc((abs(M9_RX)));
        % axis xy;
        % xlabel('OFDM sym');
        % ylabel('Subcar');
        % title('8 po sobe jdoucích SSBlocku');
        % colormap jet;
        % colorbar;

        for OFDM_symbol = 1:32
            for subcarrier_indx = 1:240
                c = M9_RX(subcarrier_indx,OFDM_symbol);
                c = c/max(max(abs(M9_RX)));

                distances = abs(c - M9_TX(subcarrier_indx,OFDM_symbol));
                sum_distances = sum(distances);

                if best_distances(subcarrier_indx,OFDM_symbol) > sum_distances
                    M9_RX_corrected(subcarrier_indx,OFDM_symbol) = c;
                    best_distances(subcarrier_indx,OFDM_symbol) = sum_distances;
                end
            end
        end

        % figure(14);
        % plot(M9_RX_corrected,'bo')
        % hold on
        % plot(M9_TX,'rx')
        % grid on
        % hold off

        sig_rx = data;
    end

    a_SSB=lscov(M9_RX_corrected,M9_TX);
    M9_RX_corrected = M9_RX_corrected*a_SSB;

    figure(15)
    %subplot(7,10,i)
    imagesc((abs(M9_RX_corrected)));
    axis xy;
    set(gca, 'FontSize', 30);
    xlabel('OFDM symbol','FontSize',35);
    ylabel('Subcarrier frequency','FontSize',35);
    title('8 Synchronization Signal Blocks','FontSize',35);
    colormap jet;
    clim([0 1])
    colorbar;

    %% Frequency spectrum
    fftSignal = fft(sig_rx);
    fftSignal = fftshift(fftSignal);
    f = fs/2*linspace(-1,1,length(fftSignal));

    % figure(16)
    % plot(f,abs(fftSignal));
    % grid on
    % title('Frequency spectrum');
    % xlabel('Frequency (Hz)');
    % ylabel('magnitude');

    %% Average power
    RSRP = abs(M9_RX).^2;
    RSRP = mean(mean(RSRP,2));
    RSRP_dBm = 10*log10(RSRP/1e-3);

    TX_direction(i) = i;
    Power_RSRP_dBm(i) = RSRP_dBm;

    %% BER and SNR
    dataMIN = abs(data(:)); % Ensure column vector

    % Determine threshold: mid-point between max and min value
    threshold = ((max(dataMIN) + min(dataMIN))*2)/3;

    % Logical masks
    high_mask = dataMIN > threshold;
    low_mask = dataMIN <= threshold;

    % Find indices of long enough high and low segments
    min_points = 10000;

    % Get continuous segments
    low_segments = get_segments(low_mask, min_points);

    % Collect samples from these segments
    low_values = collect_segment_values(dataMIN, low_segments);

    % Calculate averages
    avgMin = mean(low_values);
    avgMax = max(abs(data))*0.9;

    BER(i,1) = 20*log10(avgMax/avgMin); % SNR
    BER(i,2) = BER_fnc(M9_TX, M9_RX_corrected); % BER
end

TX_direction = TX_direction';
RP_RSRP_dBm = Power_RSRP_dBm';
T = table(TX_direction, RP_RSRP_dBm);

BER = sortrows(BER, 1);

figure(18)
semilogy(BER(:,1),BER(:,2))
grid on
xlabel('SNR [dB]');
ylabel('SER [-]');
title('Symbol error rate');
ylim([1e-5 1])

[M, I] = max(RP_RSRP_dBm);
fprintf('For the receiver position is the best beam %.1f with the power of SSBlock %.3f dBm.\r\n', smer(I), M);