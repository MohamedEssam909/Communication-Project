clc;
clear;
close all;
pkg load signal;
fs=1000; %sampling freq
ts=1/fs;
T=6.4; %simulation time
fc=100; %carrier frequency
t=0:ts:T-ts;
sym_dur=0.1;
df=1/T;
%N=fix(T/ts); %Number of time samples
sym_len=fix(fs*sym_dur);
no_sym= T/sym_dur; %64
s=randi([0 1],1,no_sym); %64 random bits 1 or 0

bin_sig=[];
for i=1:no_sym
  bin_sig=[bin_sig,repmat(s(i),1,sym_len)];
end


if(rem(sym_len,2)==0) %% Even
  f = - (0.5*fs) : df : (0.5*fs-df) ; %% Frequency vector if x/f is even
else %% Odd
  f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ; %% Frequency vector if x/f is odd
end



subplot(611)
plot(t,bin_sig)
%ylim([-0.1,1,1.5])
xlabel('Time(s)');
ylabel('amplitude');
title('Random Binary Signal');


%modulation
ask= bin_sig.*sin(2*pi*fc.*t);
subplot(612)
plot(t,ask)
ylim([0, max(ask)]);
hold on
plot(t,bin_sig)
xlabel('Time(s)');
ylabel('amplitude');
title('ASK Sigal');
%axis tight

%demodulation using Coherent detector
demod_sig=ask.*sin(2*pi*fc.*t);
demod_sig2=ask.*sin(2*pi*fc.*t+pi/6);
demod_sig3=ask.*sin(2*pi*fc.*t+pi/3);
demod_sig4=ask.*sin(2*pi*fc.*t+pi/2);
%LPF
[b,a]= butter(6,fc/fs/2); %normalization of frequency carrier/NYQ

filtered_sig_withoutPhaseAdd= filter(b,a,demod_sig);
subplot(613)
plot(t,filtered_sig_withoutPhaseAdd)
xlabel('Time(s)');
ylabel('amplitude');
title('Filtered Sigal_0');

filtered_sig_Phase30= filter(b,a,demod_sig2);
subplot(614)
plot(t,filtered_sig_Phase30)
xlabel('Time(s)');
ylabel('amplitude');
title('Filtered Sigal_30');


filtered_sig_withPhase60= filter(b,a,demod_sig3);
subplot(615)
plot(t,filtered_sig_withPhase60)
xlabel('Time(s)');
ylabel('amplitude');
title('Filtered Sigal_60');

filtered_sig_withPhase90= filter(b,a,demod_sig4);
subplot(616)
plot(t,filtered_sig_withPhase90)
ylim([0, max(filtered_sig_withPhase90)])
xlabel('Time(s)');
ylabel('amplitude');
title('Filtered Sigal_90');

%frequency domain before modulation

bin_sig_frequeny = fftshift(fft(bin_sig)) /sym_len;% 1/Number of samples for periodic signal
figure(2)
subplot(311)
plot(f,abs(bin_sig_frequeny))
xlabel('Frequency (Hz)')
ylabel('|bin_sig(f)|')
box off

%frequency domain after modulation
ASK= fftshift(fft(ask)) /sym_len;% 1/Number of samples for periodic signal
figure(2)
subplot(312)
plot(f,abs(ASK))
xlabel('Frequency (Hz)')
ylabel('|ask(f)|')
box off


%frequency domain after demodulation
FILTERED_SIG= fftshift(fft(filtered_sig_withoutPhaseAdd)) /sym_len;% 1/Number of samples for periodic signal
figure(2)
subplot(313)
plot(f,abs(FILTERED_SIG))
xlabel('Frequency (Hz)')
ylabel('|filtered_sig_withoutPhaseAdd(f)|')
box off


%frequency domain after demodulation with angles

FILTERED_SIG_30= fftshift(fft(filtered_sig_Phase30)) /sym_len;
figure(3);
subplot(311);
plot(f,abs(FILTERED_SIG_30));
xlabel('Frequency (Hz)');
ylabel('|filtered_sig_withphase_30(f)|');
box off;

FILTERED_SIG_60= fftshift(fft(filtered_sig_withPhase60)) /sym_len;
figure(3);
subplot(312);
plot(f,abs(FILTERED_SIG_60));
xlabel('Frequency (Hz)');
ylabel('|filtered_sig_withphase_60(f)|')
box off

FILTERED_SIG_90= fftshift(fft(filtered_sig_withPhase90)) /sym_len;
figure(3)
subplot(313)
plot(f,abs(FILTERED_SIG_90))
xlabel('Frequency (Hz)')
ylabel('|filtered_sig_withphase_90(f)|')
box off
