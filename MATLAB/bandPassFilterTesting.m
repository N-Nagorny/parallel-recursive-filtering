clear all;
close all;
M_min = 16; M_max = 15;
H=zeros(256, 1);
u_0 = 50;
H(129 - u_0:129 + u_0) = 1;
ax=-128:127;
figure; plot(ax, H);
h = fftshift(ifft(fftshift(H))); %ÔÐÒ ÁÈÕ-ôèëüòðà

figure; plot(ax, real(h), ax, imag(h));
if abs(max(max(imag(h)))) > 1e-10
    return;
end
h = real(h);
figure; plot(ax, h);
[h_max, x_max] = max(abs(h));
h_ideal = [0;h(x_max-M_min+1:x_max+M_max)]
figure; plot(-M_min:M_max, h_ideal);
h = zeros(256, 1);
h(x_max-M_min:x_max+M_max) = h_ideal;
figure; plot(h);
H_ideal = fftshift(fft(fftshift(h)));
ax = -128:127; ax = ax * (pi / 128);
figure; plot(ax, real(H_ideal), ax, imag(H_ideal), ax, H);


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
g(:,:,1) = conv2(h_ideal,h_ideal,f(:,:,1),'same');
g(:,:,2) = conv2(h_ideal,h_ideal,f(:,:,2),'same');
g(:,:,3) = conv2(h_ideal,h_ideal,f(:,:,3),'same');
figure; imshow(uint8(g));