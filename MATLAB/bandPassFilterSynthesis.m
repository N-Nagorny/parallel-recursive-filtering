clear all
close all
M_min = 16; M_max = 15;
N = 32;
h = zeros(32, 32);
for k = 0:31
    for m = -16:15
    h(k+1,m+17) = cos(pi / (2*N) * k * (2 * (m + M_min) + 1)) / cos(pi / (2*N) * k ) * (heaviside(m + M_min) - heaviside(m - M_max - 1));
    end
end

B = zeros(32, 32);
for l = 0:31
    for k = 0:31
        for m=-16:15
            B(l+1,k+1) = B(l+1,k+1) + h(k+1, m+17)*h(l+1, m+17);
        end
    end
end

H_inf=zeros(256, 1);
u_0 = 50;
H_inf(129 - u_0:129 + u_0) = 1;
h_inf = fftshift(ifft(fftshift(H_inf))); %ÔÐÒ ÁÈÕ-ôèëüòðà
h_inf = real(h_inf);
[h_max, x_max] = max(abs(h_inf));
h_ideal = [0;h_inf(x_max-M_min+1:x_max+M_max)]

C = zeros(32, 1);
for k = 0:31
    for m = -16:15
        C(k+1) = C(k+1) + h_ideal(m+17)* h(k+1, m+17);
    end
end

A = inv(B)*C;
A = A(1:14);
A = padarray(A, 18, 'post');
R=A.'*C;
h_new=zeros(32,1);
for m = -16:15
    for k = 0:31
        h_new(m+17)=h_new(m+17)+A(k+1)*h(k+1, m+17);
    end
end
m=-16:15;
plot(m,h_new, m,h_ideal);

f = imread('image.jpg');
figure; imshow(f);
[height, width, channels] = size(f);
noise = zeros(height, width);
u = 2; v = 0.1;
for x = 1:width
    for y = 1:height;
        noise(y, x) = 250 * cos(u * x + v * y);
    end;
end;
H=ones(256,1);

f = double(f);
f(:,:,1) = f(:,:,1) + noise;
f(:,:,2) = f(:,:,2) + noise;
f(:,:,3) = f(:,:,3) + noise;
figure; imshow(uint8(f));

g = zeros(height, width, channels);
g(:,:,1) = conv2(h_new,h_new,f(:,:,1),'same');
g(:,:,2) = conv2(h_new,h_new,f(:,:,2),'same');
g(:,:,3) = conv2(h_new,h_new,f(:,:,3),'same');
figure; imshow(uint8(g));

eps2 = 0;
for m=-16:15
    eps2 = eps2 + (h_ideal(m+17)-h_new(m+17))^2;
end

eps_min2 = 0;
for m=-16:15
eps_min2 = eps_min2 + (h_ideal(m+17))^2 - R;
end

