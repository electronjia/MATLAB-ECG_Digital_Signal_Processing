This is a project work on performing digital signal processing and data analysis on ECG data that was acquired using ECG prototype circuit (Figure 6 on report) and oscilloscope data. 

I analyzed the ECG data and performed the following:
- visualization
- Fast Fourier Transform (setting sampling frequency, getting double-sided and single-sided amplitude spectrum, signal frequencies)
- Finite Impulse Response (FIR) lowpass filter (setting Nyquist and cut-off frequencies, setting filter order)
- FIR bandpass filter (setting upper and lower bound cut-off frequencies, )
- Detecting R peaks from the QRS complex of an ECG
- Computing heart rate using 3 methods

I analyzed the ECG data with the above digital signal processing and data analysis techniques for data that used the following resistance: 4.7K Ohm, 10K Ohm, 150 Ohm, 470 Ohm

The task was to analyze the raw data that was obtained from oscilloscope and to perform digital signal processing techniques to acquire the peak frequencies, detect R peaks, and calculate heart rate. We performed the above techniques and were able to calculate the heart rate using 3 different approaches. All three approaches were close to being accurate (did not perform stats on accuracy).
