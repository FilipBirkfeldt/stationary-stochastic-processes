%% Stationära stokastiska processer 
% Datorlaboration 2 
% Filip Birkfeldt & Nils Barr Zeilon 

%% 3 
clear all; 
clc; 
% fft, periodogram, pwelch, cpsd, and mscohere

data = load('/Users/filipbirkfeldt/Desktop/StatStok/all_files/unknowndata.mat')
data = getfield(data, 'data')
plot(data)
% Mean tydligt skillt från 0. 
m = mean(data); 
%data = data-m
plot(data)
%% 
X = fft(data); 
n = length(data)
Rhat = (X.*conj(X))/n
f = [0:n-1]/n 
figure; 
plot(f, Rhat)
figure; 
plot(f-0.5,fftshift(Rhat))
% Q, What do you see? Can you say anything about the signal?
% Återkommande vid f=+/- 0.2 -> harmonic function. 
% Spectral density hög vid denna frekvens. 
figure; 
plot(data)
%% 
% Zero-padding : X=fft(x,N);
N = 4096; 
X=fft(data,N);
Rhat=(X.*conj(X))/n; 
f=[0:N-1]/N;
plot(f,Rhat)
% Q, How does the use of zero-padding change your spectral estimate? The
% unknown signal is a realization of a cosine function process, what should
% the actual covariance function look like? Can you identify the frequency?

% svar: Ser tydliga andra frekvensen och dess "power".
%       Bättre visualisering av spectral density
%       covariance-function svänger myket för tau från början till slut->avtar 

%       f=0.2 
%       
%% 
rhat = ifft(Rhat)
plot([0:15],rhat(1:16))

% Q, Can you explain the shape of the estimate from theory?
% Kovariansen är som störst vid tau=0. 
% Eftersom att det är en cosine-func vi kollar på så kommer kovariansen att
% variera mellan +/- och också få en form som en cosninus. 

%% 4
clear all; 
clc; 
e = randn(500,1);
modell.A=[1 -2.39 3.35 -2.34 0.96]; 
modell.C=[1 0 1];
x = filter(modell.C, modell.A, e);
plot(x)
figure; 
[H,w]=freqz(modell.C,modell.A,2048);
R=abs(H).^2;
plot(w/2/pi,10*log10(R))
figure; 
periodogram(x,[],4096);
figure; 
periodogram(x,hanning(500),4096);

% Q, How does the periodogram and the Hanning windowed periodogram differ? 
% Hanning har mindre bias. Klarar av dippen vid f=0.25 Hanning har dock
% bredare main-loobs - mer spectral leakage. 

%% 
[Rhat,f]=periodogram(x,[],4096,1);
figure; 
plot(f,Rhat) % Linear scale
semilogy(f,Rhat) % Logarithmic scale
% Welch
K = 10; 
n = length(x); 
L = 2*n/(K+1); 
%L=2; 
figure; 
pwelch(x,hanning(L),[],4096);
% Q. What are the differences using the Welch estimator compared to 
%    the periodogram? Sketch and explain.

% Welch metod är att ta genomsnittet av periodogram
% Därför mycket jämnare och finare kurva för Welch 

%% 
format short;
e = randn(500,1);

Rhate=periodogram(e,[],4096);
Rhatew=pwelch(e,hanning(L),[],4096);
var(Rhate)/var(Rhatew)

% Q. What is the average relation between the two variances?
% Is this in concordance with theory?

% någonstans grovt vid 1000. 
% Welch metod är kör medeltal av periodogram för att minska variansen. 
% Resultatet stödjer därför teorin.  

%% Thomsom - blev riktigt bra
pmtm(x,(10-1)/2,[],4096);

%% 5 - Spectral analysis of an EEG signal 
clear all; 
clc; 
eeg = load('/Users/filipbirkfeldt/Desktop/StatStok/all_files/eegdata12.mat')
data1 = getfield(eeg, 'data1')
data2 = getfield(eeg, 'data2')
data3 = getfield(eeg, 'data3')

% Q. Which sequence is most likely to come from the 12 Hz flickering light?
% Data. 2 

%% 
clear all; 
clc; 
eeg_dax = load('/Users/filipbirkfeldt/Desktop/StatStok/all_files/eegdatax.mat') 
data1 = getfield(eeg_dax, 'data1')
data2 = getfield(eeg_dax, 'data2')
data3 = getfield(eeg_dax, 'data3')

% Q. Can you judge which frequency it has and which sequence it is intro- duced into?
% I data3, frekvensen f=16

%% 6 - Identification of stationary Gaussian processes
clear all; 
clc; 
data = load('/Users/filipbirkfeldt/Desktop/StatStok/all_files/threeprocessdata.mat') 
% y1, y2, y3 relaizations from three different zero-mean processes. 
% white Gaussian noise sequences, x1, x2, and x3, input sequences of the
% filters.
y1 = data.y1, y2 = data.y2, y3 = data.y3;
x1 = data.x1, x2 = data.x2, x3 = data.x3;
subplot(3,1,1)
plot(y1)
title('y1')
hold on; 
subplot(3,1,2)
plot(y2)
title('y2')
hold on; 
subplot(3,1,3)
plot(y3)
title('y3')
%% 
% Periodigram på varje signal, 
% sen - kör welchmetoden på varje
% utifrån detta -> ser man vilken signal som hör ihop med vilekn spectral
% density 

c1 = xcov(y1);
c2 = xcov(y2);
c3 = xcov(y3);
% zero-padding 
n=length(y1); 
N=2048; 

Y1=fft(y1,N);
Y2=fft(y2,N);
Y3=fft(y3,N);

f=[0:N-1]/N;
Rhat1=(Y1.*conj(Y1))/n; 
Rhat2=(Y2.*conj(Y2))/n; 
Rhat3=(Y3.*conj(Y3))/n; 

figure; 
subplot(3,1,1)
plot(f,Rhat1)
title('Rhat1')
hold on; 
subplot(3,1,2)
plot(f,Rhat2)
title('Rhat2')
hold on; 
subplot(3,1,3)
plot(f,Rhat3)
title('Rhat3')

%Welch 
% L = 2*n/(K+1); 
figure; 
K=10; 
L = (2*length(y1))/(K+1)
subplot(3,1,1) 
pwelch(Y1,hanning(L),[],4096);
hold on;
subplot(3,1,2) 
pwelch(Y2,hanning(L),[],4096);
hold on;
subplot(3,1,3) 
pwelch(Y3,hanning(L),[],4096);
hold on;

% Q. Compare your periodograms and Welch spectra and the true spectral
% densities of Figure 2. Identify which sequence that belongs to each of the
% spectral densities, A, B and C. Sketch and explain.

% PASS

%% 6 
figure;
mscohere(x1,y1,hanning(L),[],4096);
figure; 
mscohere(x3,y3,hanning(L),[],4096);
%% vänder på det
figure;
mscohere(x1,y3,hanning(L),[],4096);
figure; 
mscohere(x3,y1,hanning(L),[],4096); 

% Q. What is the coherence spectral density according to theory? Relating
% to this knowledge, which of your coherence estimates seems to be correct?
% Were the names of the input sequences shifted?

% Coherence spectrum equal to 1 means the amplitudes in X′(t) are directly 
% proportional to those in X(t).
% Så det borde vara 1? alltså har dem skiftats? 






