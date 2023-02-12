phase_shift_1 = 0;%2*pi*90/360;
phase_shift_2 = 0;%2*pi*10/360;
phase_shift_3 = 0;%2*pi*30/360;
freq_shift_1 = 0;%10*2*pi;


[y1_old,FS1]=audioread("2.mp3");
[y2_old,FS2]=audioread("4.m4a");
[y3_old,FS3]=audioread("1.mp3");

y1_old = y1_old(:,1) + y1_old(:,2);
y2_old = y2_old(:,1) + y2_old(:,2);
y3_old = y3_old(:,1) + y3_old(:,2);

FS_all = 220000;
[P, Q] = rat(FS_all/FS1);
y1 = resample(y1_old,P,Q);

[P, Q] = rat(FS_all/FS2);
y2 = resample(y2_old,P,Q);

[P, Q] = rat(FS_all/FS3);
y3 = resample(y3_old,P,Q);

max_len = max(length(y1), max(length(y2), length(y3)));
t_all = linspace(0,max_len/FS_all,max_len);  

y1 = [y1;zeros(max_len-length(y1), 1)];
y2 = [y2;zeros(max_len-length(y2), 1)];
y3 = [y3;zeros(max_len-length(y3), 1)];

wc1 = bandwidth(y1) ; 
carrier1=cos(wc1*t_all);
modulated1=(y1.*carrier1');

wc2 = bandwidth(y2);
carrier2=cos(wc2*t_all);
modulated2=(y2.*carrier2');


carrier3=sin(wc2*t_all);
modulated3=(y3.*carrier3');


% Adding the modulated Signal
modulated = modulated1 + modulated2 + modulated3;


before_low_pass_demodulation1 = bandpass(modulated,[75000,125000],FS_all);
before_low_pass_demodulation2 = bandpass(modulated,[10,70000],FS_all);

%demodulation
carrier1_demodulation = cos((wc1-freq_shift_1)*t_all+phase_shift_1);
carrier2_demodulation = cos(wc2*t_all+phase_shift_2);
carrier3_demodulation = sin(wc2*t_all+phase_shift_3);

demodulation1 = before_low_pass_demodulation1.*carrier1_demodulation';
demodulation2 = before_low_pass_demodulation2.*carrier2_demodulation';
demodulation3 = before_low_pass_demodulation2.*carrier3_demodulation';


demodulation1_after_low_pass = lowpass(demodulation1, 25000,FS_all);
demodulated_signal1 = demodulation1_after_low_pass*2;
%soundsc(demodulated_signal1,FS_all);

demodulation2_after_low_pass = lowpass(demodulation2, 70000,FS_all);
demodulated_signal2 = demodulation2_after_low_pass*2;
%soundsc(demodulated_signal2,FS_all);

demodulation3_after_low_pass = lowpass(demodulation3, 70000,FS_all);
demodulated_signal3 = demodulation3_after_low_pass*2;
%soundsc(demodulated_signal3,FS_all);

disp(bandwidth(y1_old)/(2*pi));
disp(bandwidth(y2_old)/(2*pi));
disp(bandwidth(y3_old)/(2*pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%speech signal1
%plot in time domain
N22=length(y1);
t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(y1,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(1);
subplot(3,1,1);
plot(t22, y1)
title('speech signal 1')
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%speech signal2
%plot in time domain
N22=length(y2);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(y2,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(2);
subplot(3,1,1);
plot(t22, y2)
title('speech signal 2')
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%speech signal3
%plot in time domain
N22=length(y3);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(y3,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(3);
subplot(3,1,1);
plot(t22, y3)
title('speech signal 3')
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modulated1
%plot in time domain
N22=length(modulated1);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(modulated1,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(4);
subplot(3,1,1);
plot(t22, modulated1)
xlabel('Time(s)')
ylabel('Amplitude')
title('time domain of modulated signal number 1')
%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modulated2
%plot in time domain
N22=length(modulated2);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(modulated2,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(5);
subplot(3,1,1);
plot(t22, modulated2)
xlabel('Time(s)')
ylabel('Amplitude')
title('time domain of modulated signal number 2')
%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modulated3
%plot in time domain
N22=length(modulated3);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(modulated3,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(6);
subplot(3,1,1);
plot(t22, modulated3)
xlabel('Time(s)')
ylabel('Amplitude')
title('time domain of modulated signal number 3')
%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modulated
%plot in time domain
N=length(modulated);

t = (0:N-1)/FS_all; 
y_freq_old = fft(modulated,N);
F = ((-N/2):(N/2)-1).*(FS_all/N);
y_freq =fftshift(y_freq_old);
magnitudeY = abs(y_freq);        
phaseY = unwrap(angle(y_freq));  

figure(7);
subplot(3,1,1);
plot(t, modulated)
title('time domain of modulated signal')
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%after bandpass filter1
%plot in time domain
N22=length(before_low_pass_demodulation1);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(before_low_pass_demodulation1,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(8);
subplot(3,1,1);
plot(t22, before_low_pass_demodulation1);
title('after band pass filter1');
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%after bandpass filter2
%plot in time domain
N22=length(before_low_pass_demodulation2);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(before_low_pass_demodulation2,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(9);
subplot(3,1,1);
plot(t22, before_low_pass_demodulation2);
title('after band pass filter2');
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%after bandpass filter3
%plot in time domain
N22=length(before_low_pass_demodulation2);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(before_low_pass_demodulation2,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(10);
subplot(3,1,1);
plot(t22, before_low_pass_demodulation2);
title('after band pass filter3');
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot after bandpass1 then carry1
%plot in time domain
N22=length(demodulation1);
t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(demodulation1,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(11);
subplot(3,1,1);
plot(t22, demodulation1)
title('plot after bandpass1 then carry1')
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot after bandpass2 then carry2
%plot in time domain
N22=length(demodulation2);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(demodulation2,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(12);
subplot(3,1,1);
plot(t22, demodulation2)
title('plot after bandpass2 then carry2')
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot after bandpass1 then carry1
%plot in time domain
N22=length(demodulation3);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(demodulation3,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(13);
subplot(3,1,1);
plot(t22, demodulation3)
title('plot after bandpass3 then carry3')
xlabel('Time(s)')
ylabel('Amplitude')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%time domain of demodulated signal1 after lowpass filter and multipling by 2
%plot in time domain
N22=length(demodulated_signal1);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(demodulated_signal1,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(14);
subplot(3,1,1);
plot(t22, demodulated_signal1)
xlabel('Time(s)')
ylabel('Amplitude')
title('time domain of demodulated signal1 after lowpass filter and multipling by 2')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% time domain of demodulated signal2 after lowpass filter and multipling by 2
%plot in time domain
N22=length(demodulated_signal2);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(demodulated_signal2,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(15);
subplot(3,1,1);
plot(t22, demodulated_signal2)
xlabel('Time(s)')
ylabel('Amplitude')
title('time domain of demodulated signal2 after lowpass filter and multipling by 2')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% time domain of demodulated signal3 after lowpass filter and multipling by 2
%plot in time domain
N22=length(demodulated_signal3);

t22 = (0:N22-1)/FS_all; 
y_freq_old33 = fft(demodulated_signal3,N22);
F = ((-N22/2):(N22/2)-1).*(FS_all/N22);
y_freq22 =fftshift(y_freq_old33);
magnitudeY22 = abs(y_freq22);        
phaseY22 = unwrap(angle(y_freq22));  

figure(16);
subplot(3,1,1);
plot(t22, demodulated_signal3)
xlabel('Time(s)')
ylabel('Amplitude')
title('time domain of demodulated signal3 after lowpass filter and multipling by 2')

%plot in freq domain
subplot(3,1,2);

plot(F,magnitudeY22);
title('Magnitude response of signal');
ylabel('Magnitude');
subplot(3,1,3);
plot(F,phaseY22);
title('Phase response of signal');
xlabel('Frequency')
ylabel('radians');