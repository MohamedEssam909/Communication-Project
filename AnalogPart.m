pkg load communications


clc
clear
close all
fs=100; %sampling freq
ts= 1/fs; %sampling time
df=0.01; %resolution
T = 1/df; %simualation time
pi= 3.14159265359;

N= ceil (T/ts);  %no of samples

if(rem(N,2)==0)
f= (-fs/2: df: (fs/2)-df);
else
f= (-((fs/2)-df/2): df: (fs/2)-df/2);  %frequecy interval
end

t= -50: ts : ts*(N-5001);    %time interval

%x(t)
x= zeros (size (t));
x(t>-2 & t<-1)= t(t>-2 & t<-1)+2;
x(t>=-1 & t<=1)=1;
x(t>1 & t<2)= -t(t>1 & t<2)+2;

%Plot x(t)
figure (1)
plot(t,x)
title ('x in time domain');
xlabel ('t(sec)');
ylabel ('x(t)');


%Fourrier Transform of x by octave
X= fftshift (fft (x))*ts;

%Plot Fourrier Transform of x by octave
figure (2)
plot (f,abs(X));
title ('x in freq domain (by octave)');
xlabel ('freq(HZ)');
ylabel ('X(f)');

%Analytical Expression for fourrier transform of x

 XX = 4* (sinc (2*f)).^2 - (sinc (f)).^2;

XXX = 3*sinc (f).*sinc(3*f);

%Plot Analytical Expression for fourrier transform of x

figure (3)
plot (f,abs( XXX));
title ('x in freq domain (theoretically)');
xlabel ('freq(HZ)');
ylabel ('XX(f)');

figure (4)

plot (f, abs(X) , 'r', f, abs (XX)  , 'b');

title ('Comparison');
xlabel ('freq(HZ)');
ylabel ('X(f), XX(f)');
legend ('X(f)', 'XX(f)');


%BW of X(f)
Energy= sum (abs(X).^2 )*df;

E_accumulator =0;

Index = find(abs (f-0)== min (abs(f-0)));
for i= (Index : length (f))
    E_accumulator = E_accumulator + abs (X(i)).^2 *df;
    if ( E_accumulator >= 0.95*0.5*Energy)
        BW= f(i);
        break
    end
end

% LPF with BW= 1HZ
H= abs (f) <= 1;

figure (5)
plot (f,H);
title ('LPF');
xlabel ('freq(HZ)');
ylabel ('H(f)');

% Signal with LPF in freq domain
S= H.*X;
figure (6)
plot (f,abs (S));
title ('X with LPF');
xlabel ('freq(HZ)');
ylabel ('S(f)');

%Signal with LPF in time domain
s= ifft (ifftshift (S)/ts);
figure (7)
plot (t,s);
title ('x with LPF');
xlabel ('t(sec)');
ylabel ('s(t)');


% LPF with BW= 0.3HZ
HH= abs (f) <= 0.3;

figure (8)
plot (f,HH);
title ('LPF #2');
xlabel ('freq(HZ)');
ylabel ('H(f)');

% Signal with LPF in freq domain
SS= HH.*X;
figure (9)
plot (f,abs (SS));
title ('X with LPF #2');
xlabel ('freq(HZ)');
ylabel ('S(f)');

%Signal with LPF in time domain
ss= ifft (ifftshift (SS)/ts);
figure (10)
plot (t,abs(ss));
title ('x with LPF #2');
xlabel ('t(sec)');
ylabel ('s(t)');




%%%%%%%%%%%%%%%part 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_m=ceil(T/ts);
%t_m=0:ts:(N-1)*ts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=cos(2*pi*t);
m(t>6)=0;
m(t<0)=0;
figure(11)
plot(t,m)
title ('x in time domain');
xlabel ('t(sec)');
ylabel ('m(t)');

M=fftshift(fft(m))*ts ;
figure(12)
plot(f,abs(M))
title ('x in freq domain (by octave)');
xlabel ('freq(HZ)');
ylabel ('M(f)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M_analytical = 3 * sinc(2*3*(f - 1)) + 3 * sinc(2*3*(f + 1));
figure(13);
plot(f, abs(M_analytical));
title('Analytical Fourier Transform');
xlabel('Frequency (Hz)');
ylabel('|F(f)|');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (14)

plot (f, abs(M) , 'r', f, abs (M_analytical)  , 'b');

title ('Comparison');
xlabel ('freq(HZ)');
ylabel ('M(f), M_analytical(f)');
legend ('M(f)', 'M_analytical(f)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
power_total_m=sum(abs(M).^2)*df;
index_m= find(abs (f-0)== min (abs(f-0)));
power_acc_m=0;
for c_index_m = index_m:length(f)
  power_acc_m=df*abs(M(c_index_m)).^2+power_acc_m;
  if(power_acc_m>=0.95*0.5*power_total_m)
  BW_m= f(c_index_m);
  break
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%% part 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc=20;
fc2=fc+5.5;

%%%%%%%%%%%%%%%%%First signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Carrier 1 %%
c=cos(2*pi*fc*t);
figure(15)
plot(t,c)
xlabel('time(sec)')
ylabel('c(t)')


C=fftshift(fft(c))/N;
figure(16)
plot(f,abs(C))
xlabel('freq(hz)')
ylabel('C')


% LPF with BW= 1.5HZ
H= abs (f) <= 1.5;
% Signal with LPF in freq domain
X_lpf= H.*X;
%Signal with LPF in time domain
x_lpf= ifft (ifftshift (X_lpf)/ts);


%%DSB-SC Modulation%%
s1_DSB_SC=x_lpf.*c;
figure(17)
plot(t,s1_DSB_SC)
xlabel('time(sec)')
ylabel('s1(t)')

S1=fftshift(fft(s1_DSB_SC))/N;
figure(18)
plot(f,abs(S1))
xlabel('freq(hz)')
ylabel('S1')




%%Carrier 2%%
c2=cos(2*pi*fc2*t);
figure(19)
plot(t,c2)
xlabel('time(sec)')
ylabel('c2(t)')

C2=fftshift(fft(c2))/N;
figure(20)
plot(f,abs(C2))
xlabel('freq(hz)')
ylabel('C2')



%%SSB Modulation%%
s2=m.*c2;
figure(21)
plot(t,s2)
xlabel('time(sec)')
ylabel('s2(t)')

S2=fftshift(fft(s2))/N;
figure(22)
plot(f,abs(S2))
xlabel('freq(hz)')
ylabel('S2')



%%BandPass Filter%%
H_BPF2=zeros(size(f));
H_BPF2(f<(fc2) & f>(fc2-1.5))=1; %%total BW SO *0.5
H_BPF2(f>-(fc2) & f<-(fc2-1.5))=1; %%total BW SO *0
figure(23)
plot(f,abs(H_BPF2));

S2=S2.*H_BPF2;
figure(23)
plot(f,abs(S2));
xlabel('freq(hz)');
ylabel('S2 after BPF');



s2_SSB=real(ifft(ifftshift(S2)*N));
figure(24)
plot(t,s2_SSB);
xlabel('t');
ylabel('s2 after BPF');







%%%%%%%%%%%%%%%%%%%%%%FDM%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%total s%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_total=s1_DSB_SC+s2_SSB
figure(25)
plot(t,s_total)
xlabel('time')
ylabel('s(t)')



%%%%%%%%%%%%%%total s in frequency Domain%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_total=fftshift(fft(s_total))/N;
figure(28)
plot(f,abs(S_total))
xlabel('freq(hz)')
ylabel('S(f)')






%%%%%%%%%%%%%%%%%%%%Coherent demodulator%%%%%%%%%%%%%%%%%%%%%%%%%%


%first signal%

%BPF%
H_BPF_demodulator=zeros(size(f));
H_BPF_demodulator(f>(fc-1.5) & f<(fc+1.5))=1;
H_BPF_demodulator(f<-(fc-1.5) & f>-(fc+1.5))=1;


S1_afterBDF=S_total.*H_BPF_demodulator;
figure(29)
plot(f,abs(S1_afterBDF))


s1_afterBDF=ifft(ifftshift(S1_afterBDF)*N);
figure(30)
plot(t,s1_afterBDF)


%demodulation with the same carrier%
s1_afterDemodulation=2*s1_afterBDF.*c;
figure(31)
plot(t,s1_afterDemodulation)

S1_afterDemodulation=fftshift(fft(s1_afterDemodulation))/N;
figure(32)
plot(f,abs(S1_afterDemodulation))

%LPF of coherent Demodulator%
H_S1_demodulator= abs (f) <= 1.5;
figure(33)
plot(f,H_S1_demodulator)

S1_recieved=S1_afterDemodulation.*H_S1_demodulator;
figure(34)
plot(f,abs(S1_recieved))


s1_recieved=ifft(ifftshift(S1_recieved)*N);
figure(35)
plot(t,s1_recieved)
hold on
plot(t,x) %Plotting transmitted signal
legend('Recieved message','transmitted message')







%second signal%
G_BPF_demodulator=zeros(size(f));
G_BPF_demodulator(f<(fc2) & f>(fc2-1.5))=1;
G_BPF_demodulator(f>-(fc2) & f<-(fc2-1.5))=1;

S2_afterBDF=S_total.*G_BPF_demodulator;
figure(36)
plot(f,abs(S2_afterBDF))

s2_afterBDF=ifft(ifftshift(S2_afterBDF)*N);
figure(37)
plot(t,s2_afterBDF)


%demodulation with the same carrier%
s2_afterDemoduation=2*2*s2_afterBDF.*c2;
figure(38)
plot(t,s2_afterDemoduation)

S2_afterDemoduation=fftshift(fft(s2_afterDemoduation))/N;
figure(39)
plot(f,abs(S2_afterDemoduation))


%LPF of coherent Demodulator%
G_S2_demdulator= abs (f) <= 1.5;
figure(40)
plot(f,G_S2_demdulator)



S2_recieved=S2_afterDemoduation.*G_S2_demdulator;
figure(41)
plot(f,abs(S2_recieved))

s2_recieved=ifft(ifftshift(S2_recieved)*N);
figure(42)
plot(t,s2_recieved)
hold on
plot(t,m) %Plotting transmitted signal
legend('Recieved message','transmitted message')
