close all
clear all
clc

Mtotal=4267;
Mdata=240;
Mdc=0;
Ncp=0;
Nsym = 2400; % OFDM symbolu ve 20 ms, tj SSburst v max 5 ms + doplneno nulama
fmezinos = 120e3; %case D
fs = 512e6;
delka = 0.0200015625; % delka trvani není presne 20 ms, protoze pocet subnosnych na jeden OFDM symbol je zaokrouhleny z 4266.66667 na 4267
ix = 0;
nulove_sub = 4; %zvolit
moznosti = 240/nulove_sub;
smer = 1:moznosti;
s = 0; % pro subplot srovnani amplitudy + faze

for i = 31 %1:moznosti
    disp("Nacitam soubor data_rx_" + i);
    if i > 1
        clear data sig_rx sig_tx
    end
    
    filename_rx = sprintf('data/measure_20_2_2025/LOS/data_rx_%d.mat', i); 
    %filename_rx = sprintf('data/measure_20_2_2025/NLOS_directional/data_rx_%d.mat', i);
    %filename_rx = sprintf('data/measure_20_2_2025/NLOS_omni/data_rx_%d.mat', i);
    load(filename_rx); 
    filename_tx = sprintf('data/tx_data/sig_tx_%d.mat', i); 
    load(filename_tx); 
    sig_rx = data;

    %%
    [M5_TX] = OFDM_demodulator(sig_tx, Mtotal, Mdata, Mdc, Ncp, Nsym);

    figure(1)
    imagesc((abs(M5_TX)));
    axis xy;
    xlabel('OFDM sym');
    ylabel('Subcar');
    title('SS burst, demodulated from RX sig_rx');
    colorbar;

    K = 0; %konstanta pro přeskládaní ve spektru
    a_TX = (floor(size(M5_TX, 1)/2+1)-240/2)-Mdc+K;
    b_TX = (floor(size(M5_TX, 1)/2)+240/2)+Mdc+K;
    M6_TX = fftshift(M5_TX, 1);

    figure(2)
    imagesc((abs(M6_TX)));
    axis xy;
    xlabel('OFDM sym');
    ylabel('Subcar');
    title('SS burst, demodulated from RX sig_rx and rearanged with fftshift');
    colormap jet;
    colorbar;

    M7_TX = M6_TX(a_TX:b_TX,1:600);
    M7_TX = [M7_TX(1:120,:); M7_TX(120+2*Mdc+1:end,:)];

    figure(3)
    imagesc((abs(M7_TX)));
    axis xy;
    xlabel('OFDM sym');
    ylabel('Subcar');
    colormap jet;
    colorbar;

    M9_TX = M7_TX(:,[4:11, 16:23, 32:39, 44:51]);
    M9_TX = round(M9_TX,3);

    figure(4)
    imagesc((abs(M9_TX)));
    axis xy;
    set(gca, 'FontSize', 30);
    xlabel('OFDM symbol','FontSize',35);
    ylabel('Subcarrier frequency','FontSize',35);
    title('8 Synchronization Signal Blocks','FontSize',35);
    colormap jet;
    colorbar;

    %% time
    t = 0:1/fs:delka*4-(1/fs); % s
    t = t*1000; % prevod na ms
    % zobrazeni prijateho signalu
    figure(5)
    %subplot(6,10,i)
    plot(t, abs(sig_rx))
    grid("on")
    %plot(abs(sig_rx)) % bez casu bude fungovat, pokud je CP>1
    set(gca, 'FontSize', 30);
    xlabel('t [ms]','FontSize',35);
    ylabel('Amplitude [-]','FontSize',35);
    title('SS burst in time domain','FontSize',35);
    xlim([0 80])
    %% synchronizace
    % prijmout nekolik repetic celeho signalu (cely signal je 20 ms)
    % zasynchronizovat se na začátek, vyseknout 20 ms a potom demodulovat
    [korelace, lags] = xcorr(abs(sig_rx), abs(sig_tx));

    figure(6);
    %subplot(6,10,i)
    plot(lags, korelace);

    %zde vytahnout z sig_rx zpet puvodni delku pro 20 ms
    [M, I] = max(korelace);
    k = 0;

    if lags(I)+k+length(sig_tx) > 4*length(sig_tx)
        I = I - length(sig_tx);
    end

    % sig_rx = [sig_rx(lags(I)+1+k:lags(I)+k+length(sig_tx))]; vyseknuti 20ms
    % okna
    sig_rx = [sig_rx(lags(I)+1+k:end)]; % vyseknuti synchronizovaného okna od
    % začátku synchronizace až do konce, kvůli přítomnosti CP

    %sig_rx = sig_rx/max(abs(sig_rx)); % normalizace

    %% srovnani synchronizace
    % s = s + 1;
    % figure(3);
    % subplot(12,10,s)
    % plot(abs(sig_rx)); hold on;
    % plot(abs(sig_tx));
    % %legend('RX', 'TX');
    % title('Time'); hold off;
    % s = s + 1;
    % subplot(12,10,s)
    % plot(angle(sig_rx)); hold on;
    % plot(angle(sig_tx));
    % %legend('RX phase', 'TX phase');
    % title('Phase'); hold off;

    %% demodulator
    [M5] = OFDM_demodulator(sig_rx, Mtotal, Mdata, Mdc, Ncp, Nsym);

    figure(7)
    %subplot(6,10,i)
    imagesc((abs(M5)));
    %clim([0 4]);
    axis xy;
    xlabel('OFDM sym');
    ylabel('Subcar');
    % title('SS burst, demodulated from RX sig_rx');
    colorbar;

    % % a přeskládání zpět do formy- data kolem nuly
    % s = floor(size(M5, 1)/2);
    K = 3; %konstanta pro přeskládaní ve spektru
    a = (floor(size(M5, 1)/2+1)-240/2)-Mdc+K;
    b = (floor(size(M5, 1)/2)+240/2)+Mdc+K;
    % M6 = [M5(Mdata/2+1:s,:); M5(1:Mdata/2,:); M5(end-Mdata/2+1:end,:); M5(s+1:end-Mdata/2,:)];

    %% fftshift - otocení spektra zpět, aby bylo symetrické kolem nuly
    M6 = fftshift(M5, 1);
    figure(8)
    %subplot(6,10,i)
    imagesc((abs(M6)));
    axis xy;
    xlabel('OFDM sym');
    ylabel('Subcar');
    % title('SS burst, demodulated from RX sig_rx and rearanged with fftshift');
    colormap jet;
    colorbar;
    %% dekodovani smeru
    M7 = M6(a:b,1:600); %detail oblasti, kde je SSBurst
    M7 = [M7(1:120,:); M7(120+2*Mdc+1:end,:)];
    %PBCH je na pozici 5,9,17,21, 33, 37, 45, 49...
    %projit cely tento OFDM symbol, kde bude na techno 24 pozicich nejnizsi celkova
    %energie, tak urcim jako danou pozici

    %vylepseni...vytipovat po sobe jdouci PBCH OFDM symboly a
    %prumerovat sumy na tech pozicich, mohlo by zpresnit vysledky
    M8 = M7(:,[5,9,17,21, 33, 37, 45, 49]);

    figure(9)
    %subplot(6,10,i)
    imagesc((abs(M7)));
    axis xy;
    xlabel('OFDM sym');
    ylabel('Subcar');
    % title('prvni PBCH Z osmi po sobe jdoucich SSBlock');
    % colormap jet;
    colorbar;
    %%
    n = 0;
    k = 0;
    HS = zeros(moznosti,8); % HS...hledani smeru, suma hodnot
    for j = [5 9 17 21 33 37 45 49]
        k = k + 1;
        for x = 1:nulove_sub:240
            n = n + 1;
            HS(n,k) = sum(20*log10(abs(M8(x:x+nulove_sub-1))));
        end
        n = 0;
    end
    HS = mean(HS,2); % prumer z 8 PBCH
    %
    % figure(7)
    % subplot(6,10,i)
    % bar(HS);
    % xlabel('pozice');
    % ylabel('vykon');
    % % title('Intenzita subnosnych v PBCH');

    [M, I] = min(HS);% I odpovida pozici a pozice odpovida smeru
    %% Reference Signal Received Power (RSRP)
    % naskladam za sebe 8 SSBlock symbolu a spocitam pro ne RSRP
    M9 = M7(:,[4:11, 16:23, 32:39, 44:51]);% M7(:,5);
    figure(10)
    %subplot(6,10,i)
    imagesc((abs(M9)));
    axis xy;
    set(gca, 'FontSize', 30);
    xlabel('OFDM symbol','FontSize',35);
    ylabel('Subcarrier frequency','FontSize',35);
    title('8 Synchronization Signal Blocks','FontSize',35);
    colormap jet;
    colorbar;

    RSRP = abs(M9).^2;
    RSRP = mean(mean(RSRP,2));
    RSRP_dBm = 10*log10(RSRP/1e-3);


    %%
    % vyhodnoceni odhadu smeru
    ix = ix + 1;
    if smer(ix) == I
        eval = 1;
    else
        eval = 0;
    end

    TX_direction(ix) = smer(ix);
    Estimated_direction(ix) = I;
    Evaluation(ix) =eval; %logical([eval]);
    Power_RSRP_dBm(ix) = RSRP_dBm;
end
%%
Evaluation = logical(Evaluation);
TX_direction = TX_direction';
Estimated_direction = Estimated_direction';
Evaluation = Evaluation'; %logical([eval]);
Power_RSRP_dBm = Power_RSRP_dBm';
RP_RSRP_dBm = Power_RSRP_dBm;
%RP_RSRP_dBm = RP_RSRP_dBm/max(RP_RSRP_dBm);
T = table(TX_direction,Estimated_direction, RP_RSRP_dBm, Evaluation)

%ze vsech Power_RSRP_dBm vyhodnotit smer s nejvyssim vykonem

[M, I] = max(RP_RSRP_dBm);

fprintf('For the receiver position is the best beam %.1f with the power of SSBlock %.3f dBm.\r\n', smer(I), M);