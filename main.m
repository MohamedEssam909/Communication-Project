
stream = randi([0 , 1],1,64);
Tb=0.5;
[t , f , x ]=manchester(stream , Tb);
N=length(x);
figure(1)
plot(t,x)
ylabel('pulse')
xlabel('t')


spectrum= fftshift(abs(fft(x)));
figure(2)
plot(abs(f),(spectrum*2)/N);
ylabel('spectrum')
xlabel('f')
ylim([0 , 1])
xlim([0 , 10])



