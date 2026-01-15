## Results: Range, Velocity, and Angle Estimation (AWR2243)

The following results are obtained from raw FMCW radar data collected using the Texas Instruments AWR2243 mmWave radar sensor. The implemented signal processing pipeline estimates the **range**, **velocity**, and **angle of arrival (AoA)** of detected targets.

---

### Detected Targets Summary

<img width="412" height="183" alt="Detected Targets Table" src="https://github.com/user-attachments/assets/26957fb4-9e3e-4908-aec1-5a657a0e08f9" />

The table above summarizes the radar detections:

- Total number of detected targets: **7**
- For each target, the estimated:
  - **Range (m)**
  - **Radial velocity (m/s)**
  - **Angle of arrival (degrees)**

This quantitative output is derived after peak detection on the range–Doppler and angle spectra.

---

### Range–Velocity Estimation

<img width="411" height="258" alt="Range Velocity Estimation" src="https://github.com/user-attachments/assets/98f889c2-31ff-4101-bba2-1d984b8f891e" />

This figure shows the **joint Range–Velocity estimation** of detected targets.

- X-axis: **Velocity (m/s)**
- Y-axis: **Range (m)**
- Z-axis: **Signal magnitude / power**
- Peaks correspond to detected targets in range–Doppler space

This result is obtained by applying:
- FFT along the fast-time dimension for range estimation
- Doppler FFT along the slow-time (chirp) dimension for velocity estimation

---

### Range–Angle (Azimuth) Estimation

<img width="411" height="252" alt="Range Angle Estimation" src="https://github.com/user-attachments/assets/d5b2823f-4d89-4c63-81ea-50a2c42e162e" />

This polar plot represents the **Range–Azimuth (Angle of Arrival) estimation**.

- Radius: **Range (m)**
- Angle: **Azimuth angle (degrees)**
- Detected peaks indicate the angular position of targets

Angle estimation is performed using spatial FFT / beamforming across the multiple receive antennas of the AWR2243 sensor.

---

### Summary

- FMCW radar signal processing successfully implemented  
- Accurate estimation of **range**, **velocity**, and **angle**  
- Validated using real-world data from the AWR2243 sensor  

These results demonstrate the effectiveness of the processing pipeline for applications such as automotive radar, target localization, and autonomous sensing.
