%% Stationära stokastiska processer 
% Datorlaboration 2 
% Filip Birkfeldt & Nils Barr Zeilon 

%% 3 
% A(z−1)=1+a1z−1 +···+apz−p
% C(z−1)=1+c1z−1 +···+cqz−q


a1 = 0.5; 
C = 1; 
A=[1, a1]; 

[H,w]=freqz(C, A); 
R=abs(H).^2;


% Plots the specteal density 
figure; 
plot(w/2/pi,R)

%Calculate covariance-function 
H=freqz(C,A,512,"whole");
Rd=abs(H).^2;
r=ifft(Rd);
figure;
stem([0:49],r(1:50))

% Q.1 Draw rough sketches of the covariance functions for a1 < 0 and a1 > 0. How do they differ?
% Ett posivt värde på a1 ger skiftande +/- värden på r(tau) 
% Ett negativt värde på a1 ger bara positiva värden på r(tau)

% Q.2 Sketch the spectral densities. At what frequency does the spectral
% density reaches the maximum value for a1 < 0 and a1 > 0, respectively?
% a1 = 0.5 ger största värde på f=0.5 
% a2 = -0.5 ger största värde på f=0

% Fråga varför det blir som det blir. 

%%

m=0;
sigma=1; 
n=400; 
e=normrnd(m,sigma,1,n); 
%%
figure; 
plot(e); 
figure; 

a1 = 0.5; 
C = 1; 
A=[1, a1]; 

x = filter(C, A, e);
plot(x)

% Q. Sketch the view of the realizations. Express Xt in Xt−1 and et for the
% two cases and explain the differences in how the realizations evolve with
% time.

% a>0 -> trycker ihop punkterna och gör punkterna mindre. 
% störtsta värdet till  ~3.444 till ~3.14

% a<0 dra upp avstånden mellan punkterna. 

% Xt = et +(-)a1*Xt-1

%% 
A = [1, -1, 0.5]
C = [1]
P = roots(A)
Z = roots(C)

figure; 
zplane(Z,P)
figure; 
[H,w]=freqz(C, A); 
R=abs(H).^2;
% Plots the specteal density 
plot(w/2/pi,R)

%Q. Which parts of the spectral density have high power and which has low 
% power respectively in reference to the location of poles and zeros?

% Låga frekveser -> -0.25 o lägre ger low power
% högre frekvenser än 0.25 ger low power
% frekvenser emellan -0.25-0.25 ger hög power med peak f=0.1 

%% 4 Study of ARMA-processes using armagui

% Q. How do the covariance function and spectral density change when you
% increase the angle? How do they change when you increase the distance
% from the unit circle?

% Desto större frekvens desto större svängningar 
% Increase distance -> tar längre tid för kovariansen att gå mot noll. 

% 
% Q. Sketch the spectral density and covariance function. Consider the
% model for an AR-process, can you guess what the ak coefficients are?
% (Hint: do the spectral density and covariance function look familiar?)

% Högre avstånd -> högre & tydligare peak. 
% ak-coeff motsvar på vilken frekvens man har peaks på. 

% börjar på värden mellan -2,2 sen blir superstor/små värden. 
% eftersom att värdet beror på det föreliggande värdet och multipliceras 
% med värden med större än 1 -> växer exp. 

%% 4.2
% Q. For which values of τ do the covariance functions of the different 
% processes above become zero? How is this value of τ related to the order?

% MA2 - tau=3 
% MA3 - tau=5
% MA4 - tau=7
% MA5 - tau=9 
% dubbla MA -1 ger det tau som gör kovariansen till 0. 

% Q. Sketch two spectral densities, one with complex zeros close to the unit
% circle and a second with complex zeros close to the origin. Explain how
% the locations of the zeros affect the spectral density.

% Nära enhetscirkeln -> går tydlig mot noll i R(f) 
% Nära origo -> ej nära noll, väldigt liten sänka i R(f)

%% 5 Speech modeling in cell phones

% Q. How much storage is saved by sending the AR-parameters 
% instead of the audio itself?

% 87% plats är sparat 142336 bytes 











