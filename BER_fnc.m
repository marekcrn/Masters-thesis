function [BER] = BER_fnc(data_TX, data_RX)
%% Splitting resource block into separate parts (PBCH, PSS, SSS)
PBCH1_TX = data_TX(:, [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32]);
PBCH2_TX = data_TX([1:47,193:240], [3,7,11,15,19,23,27,31]);
PBCH_TX = [PBCH1_TX(:); PBCH2_TX(:)];

PSS_TX = data_TX(57:183,[1,5,9,13,17,21,25,29]);
PSS_TX = PSS_TX(:);

SSS_TX = data_TX(57:183,[3,7,11,15,19,23,27,31]);
SSS_TX = SSS_TX(:);

PBCH1_RX = data_RX(:, [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32]);
PBCH2_RX = data_RX([1:47,193:240], [3,7,11,15,19,23,27,31]);
PBCH_RX = [PBCH1_RX(:); PBCH2_RX(:)];

PSS_RX = data_RX(57:183,[1,5,9,13,17,21,25,29]);
PSS_RX = PSS_RX(:);

SSS_RX = data_RX(57:183,[3,7,11,15,19,23,27,31]);
SSS_RX = SSS_RX(:);

%% Constellation diagrams
figure(17)
subplot(1,2,1)
hold on
grid on
scatter(real(PBCH_TX), imag(PBCH_TX), 100, 'o', 'red', "filled");
scatter(real(SSS_TX), imag(SSS_TX), 100, 'o', 'green', "filled");
scatter(real(PSS_TX), imag(PSS_TX), 100, 'o', 'blue', "filled");
xlim([-1 1])
ylim([-1 1])
set(gca, 'FontSize', 30);
legend({'PBCH', 'SSS', 'PSS'},'FontSize',35)
xlabel('I','FontSize',35);
ylabel('Q','FontSize',35);
title('IQ diagram for transmitted signal','FontSize',35);

subplot(1,2,2)
hold on
grid on
scatter(real(PBCH_RX),imag(PBCH_RX), 100,'red', "filled")
scatter(real(SSS_RX),imag(SSS_RX), 100,'green', "filled")
scatter(real(PSS_RX),imag(PSS_RX), 100,'blue', "filled")
xlim([-1 1])
ylim([-1 1])
set(gca, 'FontSize', 30);
legend({'PBCH', 'SSS', 'PSS'},'FontSize',35)
xlabel('I','FontSize',35);
ylabel('Q','FontSize',35);
title('IQ diagram for recieved signal','FontSize',35);

%% PBCH Symbols
symbols_PBCH = [complex(0.530,0.530), complex(-0.530,0.530), complex(0.530,-0.530), complex(-0.530,-0.530), ...
                complex(0.707,0.707), complex(-0.707,0.707), complex(0.707,-0.707), complex(-0.707,-0.707)]; 

for k = 1:length(PBCH_RX)
    PBCH_RX(k) = symbols_PBCH(findClosestSymbol(PBCH_RX(k), symbols_PBCH));
end

%% PSS Symbols
symbols_PSS = [0.25, -0.25]; 

for k = 1:length(PSS_RX)
    PSS_RX(k) = symbols_PSS(findClosestSymbol(PSS_RX(k), symbols_PSS));
end

%% SSS Symbols
symbols_SSS = [0.5, -0.5]; 

for k = 1:length(SSS_RX)
    SSS_RX(k) = symbols_SSS(findClosestSymbol(SSS_RX(k), symbols_SSS));
end

%% Bit Error Rates
BER_PBCH = sum(PBCH_TX ~= PBCH_RX);
BER_PSS = sum(PSS_TX ~= PSS_RX);
BER_SSS = sum(SSS_TX ~= SSS_RX);

BER = (BER_PBCH + BER_PSS + BER_SSS) / (length(PBCH_TX) + length(PSS_TX) + length(SSS_TX));
end

%% Helper Function: Find Closest Symbol
function idx = findClosestSymbol(value, symbolArray)
[~, idx] = min(abs(symbolArray - value));
end
