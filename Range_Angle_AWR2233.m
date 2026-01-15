function HardWare
 chirpParams.slope = 60.012e12; % 60 MHz per µs
 chirpParams.T = 30e-6;
 chirpParams.f0 = 77e9;
 chirpParams.repeat_duration = 160e-6;
 chirpParams.B = chirpParams.slope * chirpParams.T;
 
 radarparams.fs = 10e6;
 radarparams.Ts = 1 / radarparams.fs;
 radarparams.N = 256;
 radarparams.L = 128;
 radarparams.P = 4;
 
 N = radarparams.N; 
 L = radarparams.L; 
 P = radarparams.P;
 NumFrames = 8;
 
 X = parseMultiChannelADCData('adc_data.bin', N, L, P, NumFrames);
 samplesPerFrame = N * L;
 X = X(1:samplesPerFrame, :);
 AdcDataTensor = zeros(N, L, P);
 for i = 1:P
     tensor = reshape(X(:, i), [N, L]);
     AdcDataTensor(:, :, i) = tensor;
 end
 DetectVisualizeTargets(radarparams, chirpParams, AdcDataTensor);
end


function [adc_mc_cmplx] = parseMultiChannelADCData(filename, N, L, P, NumFrames)
FID = fopen(filename, 'r');
data = fread(FID, 'int16');
fclose(FID);

adc_mc_cmplx = zeros(N * L * NumFrames, P);
for ch = 1:P
    ch_offset = 2 * (ch - 1) + 1;
    for sampIdx = 1:N * L * NumFrames
        im = data(2 * P * (sampIdx - 1) + ch_offset);
        re = data(2 * P * (sampIdx - 1) + ch_offset + 1);
        adc_mc_cmplx(sampIdx, ch) = complex(re, im);
    end
end
end

function DetectVisualizeTargets(rp, cp, x_adc) 
% 
======================================================================== 
% Function: DetectVisualizeTargets 
% ------------------------------------------------------------------------ 
% Description: 
%   
Performs Range–Velocity–Angle estimation using 3D FFT processing 
%   
from raw multi-channel ADC data. 
%   
1. Range FFT (fast time) 
%   
2. Doppler FFT (slow time) 
%   
3. Angle FFT (spatial) 
%   
Includes CFAR detection and Range–Azimuth visualization. 
% 
======================================================================== 
 
c = 3e8;  % Speed of light (m/s) 
N = rp.N; L = rp.L; P = rp.P; 
DFT_3 = zeros(N, L, P); 
 
%% ------------------- RANGE–DOPPLER PROCESSING -------------------------- 
% Perform 2D FFT for each antenna channel (Range & Doppler) 
for p = 1:P 
   mat = x_adc(:, :, p);          % Extract antenna channel data 
   rowfft = fft(mat, N, 1);       % FFT along range (fast time) 
   DFT2D = fft(rowfft, L, 2);     % FFT along Doppler (slow time) 
   DFT2D = fftshift(DFT2D, 2);    % Shift Doppler zero frequency to 
center 
   DFT_3(:, :, p) = DFT2D;        % Store result for current antenna 
end 
 
%% ------------------- ANTENNA COMBINATION ------------------------------- 
% Non-coherent accumulation (sum of magnitudes) across antennas 
DFT_accum = zeros(N, L); 
for i = 1:P 
   DFT_accum = DFT_accum + abs(DFT_3(:, :, i)); 
end 
 
% Power map (used for CFAR detection) 
DFT_power = sum(abs(DFT_3).^2, 3); 
 
%% ------------------- RANGE AND VELOCITY AXES --------------------------- 
raxis = 0 : (c * cp.T) / (2 * cp.B * rp.N * rp.Ts) : ... 
       (c * cp.T) / (2 * cp.B * rp.N * rp.Ts) * (rp.N - 1); 
vaxis = (c / (2 * rp.L * cp.f0 * cp.repeat_duration)) * ... 
       (-rp.L / 2 : rp.L / 2 - 1); 
[Vaxis, Raxis] = meshgrid(vaxis, raxis); 
 
%% ------------------- RANGE–VELOCITY MAP ------------------------------- 
figure; 
surf(Vaxis, Raxis, abs(DFT_accum), 'EdgeColor', 'none'); 
colorbar; 
xlabel('Velocity (m/s)'); 
ylabel('Range (m)'); 
title('Range and Velocity Estimation'); 
 
%% ------------------- CFAR DETECTION ----------------------------------- 
DFT2 = DFT_power; 
[Nr, Lr] = size(DFT2); 
detections = zeros(Nr, Lr); 
 
threshold_dB = 35; 
NoiseBandWidth = 1; 
GuardBandWidth = 1; 
threshold = 10^(threshold_dB / 10); 
halfWin = NoiseBandWidth + GuardBandWidth; 
 
for r = 1 + halfWin : Nr - halfWin 
   for cInd = 1 + halfWin : Lr - halfWin 
       window = abs(DFT2(r - halfWin : r + halfWin, ... 
           cInd - halfWin : cInd + halfWin)).^2; 
       % Exclude guard band (immediate neighbors) 
       window(halfWin + 1 - GuardBandWidth : halfWin + 1 + 
GuardBandWidth, ... 
              halfWin + 1 - GuardBandWidth : halfWin + 1 + 
GuardBandWidth) = 0; 
       noiseAvg = mean(window(window > 0)); 
       CUT = abs(DFT2(r, cInd))^2; 
       if CUT > threshold * noiseAvg 
           detections(r, cInd) = 1; 
       end 
   end 
end 
 
%% ------------------- CFAR RESULT VISUALIZATION ------------------------- 
figure; 
surf(Vaxis, Raxis, double(detections), 'EdgeColor', 'none'); 
colorbar; 
xlabel('Velocity (m/s)'); 
ylabel('Range (m)'); 
title('CFAR Detections (Range-Velocity)'); 
view(2); 
 
%% ------------------- EXTRACT DETECTED PEAKS ---------------------------- 
peak_id = 1; 
rInd = []; 
vInd = []; 
for r = 1:size(detections, 1) 
   for cInd = 1:size(detections, 2) 
       if detections(r, cInd) == 1 
           fprintf('Peak %d: Range = %.2f m, Velocity = %.2f m/s\n', ... 
               peak_id, raxis(r), vaxis(cInd)); 
           rInd(peak_id) = r; 
           vInd(peak_id) = cInd; 
           peak_id = peak_id + 1; 
       end 
   end 
end 
 
Rest = Raxis(rInd); 
angles = []; 
 
%% ------------------- ANGLE ESTIMATION (3D FFT across antennas) --------- 
for i = 1:length(rInd) 
   angle_v = squeeze(DFT_3(rInd(i), vInd(i), :));   % Extract antenna 
dimension 
   fft_angle = fftshift(fft(angle_v, P*4));         % 4x interpolation 
for finer angle resolution 
   lambda = c / cp.f0; 
   [~, angleindex] = max(abs(fft_angle));           % Find max peak 
   figure; 
   angleAxis = [-((P*4)/2): 1 : ((P*4)/2-1)];       % FFT bins 
   angleAxis = asin((2/(4*P))*angleAxis)*(180/pi);  % Map FFT bins to 
degrees 
   angles(i) = angleAxis(angleindex); 
   plot(angleAxis, 20 * log10(abs(fft_angle))); 
   xlabel('Angle (degrees)'); 
   ylabel('Normalized Magnitude (dB)'); 
   title('Angle Spectrum (3D-DFT across Antenna Array)'); 
   grid on; 
end 
 
%% ------------------- RANGE–AZIMUTH POLAR PLOT -------------------------- 
figure; 
disp(angles); 
polarscatter(deg2rad(angles), Rest, 'filled');  % Convert to radians for 
polar plot 
title('Range–Azimuth Polar Plot'); 
set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise'); 
thetalim([-90 90]); 
 
ax = gca; 
rlim([0 max(Rest)+1]); 
rticks(0:5:max(Rest)+1); 
 
text(pi/2, ax.RLim(2) * 1.05, 'Range (m)', ... 
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 
'FontWeight', 'bold'); 
text(0, -ax.RLim(2) * 0.15, 'Angle (°)', ... 
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 
'FontWeight', 'bold'); 
 
end
