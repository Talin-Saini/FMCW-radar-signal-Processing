function [xIF_noisy] = GenerateSingleTargetMultiChirpIFSignal(rp,cp,R,v)
c = 3e8;
adc_start_time = rp.adc_start;
t = (adc_start_time:rp.Ts:(rp.N-1)*rp.Ts+adc_start_time)';
xIF_noisy = zeros(rp.N, rp.L);
for index = 1:rp.L
 phaseTx = (2*pi*cp.f0*t) + (pi*cp.slope*(t.^2));
 tau = 2*(R + (index-1)*cp.repeat_duration*v + v*t)/c;
 phaseRx = (2*pi*cp.f0*(t-tau)) + (pi*cp.slope*(t-tau).*(t-tau));
 phaseDiff = phaseTx - phaseRx;
 xIF_noisy(:,index) = complex(cos(phaseDiff), sin(phaseDiff));
end
end
function[xIF_noisy] = GenerateMultiTargetIFSignal(rp,cp,Rvec,Vvec)
c = 3e8;
N = rp.N;
L = rp.L;
XIF_noisy = zeros(N, L);
for i = 1:length(Rvec)
 phase = 2*pi*rand();
 attenuation_factor = complex(cos(phase),sin(phase));
 XIF_noisy = XIF_noisy +
attenuation_factor*GenerateSingleTargetMultiChirpIFSignal(rp, cp,
Rvec(i), Vvec(i));
end
xIF = XIF_noisy;
noisevar = (10^(rp.SNR/10))^-1;
real_noise = sqrt(noisevar/2)*randn(N,L);
compl_noise = sqrt(noisevar/2)*randn(N,L);
w = complex(real_noise, compl_noise);
xIF_noisy = xIF + w;
end
function[DFT2] = CalcVisualize2dDFT(xIF_noisy, rp, cp)
DFT_1 = zeros(rp.N,rp.L);
c = 3e8;
for index = 1:rp.L
 %dft_1 = dft(xIF_noisy(:,index));
 DFT_1(:,index) = fft(xIF_noisy(:,index));
end
DFT_2 = zeros(rp.N,rp.L);
for index = 1:rp.N
 %dft_2 = dft(DFT(index,:));
 DFT_2(index,:) = fftshift(fft((DFT_1(index,:))));
end
%DFT_2 = fftshift(DFT_2);
DFT2 = DFT_2;
% Define velocity and range axis for plotting
raxis =
[0:(c*cp.T)/(2*cp.B*rp.N*rp.Ts):(c*cp.T)/(2*cp.B*rp.N*rp.Ts)*(rp.N-1)];
vaxis = [c/(2*rp.L*cp.f0*cp.repeat_duration)*(-rp.L/2):
c/(2*rp.L*cp.f0*cp.repeat_duration):c/(2*rp.L*cp.f0*cp.repeat_duration)*
(rp.L/2-1)];
[Vaxis, Raxis] = meshgrid(vaxis, raxis);
figure;
surf(Vaxis, Raxis, (abs(DFT_2)), 'EdgeColor','none');
colorbar;
% Display the estimated range and velocity
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Range and Velocity Estimation');
end
Q2)
function detections = CFAR_detector(DFT2, threshold_dB,
NoiseBandWidth, GuardBandWidth, rp, cp)
c = 3e8; %
[N,L] = size(DFT2);
detections = zeros(N,L);
threshold = 10^(threshold_dB/10);
halfWin = NoiseBandWidth + GuardBandWidth;
for r = 1+halfWin : N-halfWin
for cInd = 1+halfWin : L-halfWin
% Local window
window = abs(DFT2(r-halfWin:r+halfWin,
cInd-halfWin:cInd+halfWin)).^2;
window(halfWin+1-GuardBandWidth : halfWin+1+GuardBandWidth,
...
halfWin+1-GuardBandWidth : halfWin+1+GuardBandWidth)
= 0;
noiseAvg = mean(window(window > 0));
CUT = abs(DFT2(r,cInd))^2;
if CUT > threshold * noiseAvg
detections(r,cInd) = 1;
end
end
end
raxis =
[0:(c*cp.T)/(2*cp.B*rp.N*rp.Ts):(c*cp.T)/(2*cp.B*rp.N*rp.Ts)*(rp.N1)];
vaxis = [c/(2*rp.L*cp.f0*cp.repeat_duration)*(-rp.L/2):
c/(2*rp.L*cp.f0*cp.repeat_duration):c/(2*rp.L*cp.f0*cp.repeat_durat
ion)*(rp.L/2-1)];
[Vaxis, Raxis] = meshgrid(vaxis, raxis);
% Plot detections
figure;
surf(Vaxis, Raxis, double(detections), 'EdgeColor','none');
colorbar;
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('CFAR Detections (Range-Velocity)');
peak_id = 1;
for r = 1:size(detections,1)
for cInd = 1:size(detections,2)
if detections(r,cInd) == 1
fprintf('Peak %d: Range = %.2f m, Velocity = %.2f
m/s\n', ...
peak_id, raxis(r), vaxis(cInd));
peak_id = peak_id + 1;
end
end
end
end

d
B
