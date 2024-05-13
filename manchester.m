function [t, f, x] =manchester(bits, Tb)

  T=length(bits)*Tb;   %total time
  n=200; %number of sampels per bit
  N=n*length(bits);  %no of total samples in whole stream duration
  dt=T/N;  %time step
  fs=1/dt; %sampeling frequency
  df=1/T;   %frequency step
  t=0:dt:T-dt;  %=N points
  x=zeros(1,length(t));

  if(rem(N,2)==0) %% Even
  f = - (0.5*fs) : df : (0.5*fs-df) ; %% Frequency vector if x/f is even
  else %% Odd
  f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ; %% Frequency vector if x/f is odd
  end

for j=0:length(bits)-1
  %one representation
if bits(j+1)==1
   x(j*n+1:(j*n)+(n/2))=1;         %1/2 bit duration positive value
   x((j*n)+(n/2)+1:(j*n)+n)=-1;    %1/2 bit duration -ve value

  %zero representation
else
   x((j*n)+1:(j*n)+(n/2))=-1;       %1/2 bit duration -ve value
   x((j*n)+(n/2)+1:(j*n)+n)=1;      %1/2 bit duration positive value

end
end
end



