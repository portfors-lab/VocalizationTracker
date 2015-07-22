y = randn(1,2^20);
pinkNoise = Pinkify(y);
sampleRate = 384610/2;
welch(y,sampleRate,0.001)
pinkNoise = Pinkify(y);
welch(pinkNoise,sampleRate,0.001)
set(gca,'yScale','log')
set(gca,'xScale','log')

