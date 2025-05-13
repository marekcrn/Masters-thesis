function [M3] = OFDM_demodulator(sRx, Mtotal, Mdata, Mdc, Ncp, Nsym)

%Mzeros=Mtotal-Mdata-2*Mdc;
Lframe=Nsym*(Mtotal+Ncp);
frame=sRx(1:Lframe);
M1=reshape(frame,Mtotal+Ncp,Nsym);
M2=M1(Ncp+1:end,:);
M3=fft(M2);
end

