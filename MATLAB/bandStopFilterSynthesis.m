close all;
clear all;
Mm=16;Mp=15;
H=ones(256,1);
u0=20; %частота среза фильтра
H(129-u0:129+u0)=0;
ax=-128:127;
figure; plot(ax, H);

h=fftshift(ifft(fftshift(H))); 
if (abs(max(max(imag(h))))>1e-10) return; end
h=real(h); [hmax,xmax]=max(abs(h));

ax=ax*(pi/128);
figure; plot (ax,h);

hidl=[0;h(xmax-Mm+1:xmax+Mp)]; %идеальна€‘–“, которую нужно синтезировать
figure; plot (hidl);
N=Mp+Mm+1;
for j=1:N
h_k(j)={@(m) cos(pi.*(j-1).*(2.*(m+Mm)+1)./(2.*N)) ./cos(pi.*(j-1)./(2.*N)).*(heaviside(m+Mm)-heaviside(m-Mp-1))}; %базис
end;
B=zeros(N,N);
for l=1:N
hl=h_k{l};
for k=1:N
hk=h_k{k};
for m=-Mm:Mp
B(k,l)= B(k,l) + hl(m)*hk(m);
end;
end;
end;
C=zeros(N,1);
for k=1:N
hk=h_k{k};
for m=-Mm:Mp
C(k,1)= C(k,1) + hidl(m+Mm+1)*hk(m);
end;
end;

A=(B^-1)*C;
hhh=@(m) (0);
for k=1:N
if (ismember(k,[1,2,3,4,5,6]) == 0)    
hk=h_k{k};
hhh=@(m) (hhh(m)+A(k,1)*hk(m));
end;
end;

m=-Mm:Mp;
figure; plot (m,hidl(m+Mm+1),m,hhh(m));
