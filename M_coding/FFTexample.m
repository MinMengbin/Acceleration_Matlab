t = 0:.001:.25;
x = sin(2*pi*50*t) + sin(2*pi*120*t);
y = x + 2*randn(size(t));
plot(y);
title('Noisy time domain signal');