%% Stationära stokastiska processer 
% Datorlaboration 1 
% Filip Birkfeldt & Nils Barr Zeilon 

filename = '/Users/filipbirkfeldt/Desktop/StatStok/all_files/data1.mat';
import_data = importdata(filename);
data = getfield(import_data, 'data') ;
data1 = getfield(import_data, 'data1'); 

%% 3 
plot(data1.x);
medel_data = mean(data1.x);
disp(medel_data);
% - *It looks like it has about zero mean yes. However the mean of the current realisation is 0.203876184400000.*
%% 
% För 95% konfidensintervall 
sigma = std(data1.x);
kvantil = 1.96; 
I_high = medel_data + (sigma*kvantil)/sqrt(length(data))
I_low = medel_data - (sigma*kvantil)/sqrt(length(data))
% A 95% confidence interval for the mean of the realisation: [-1.700348632407131,2.108101001207130]*
%% 3.2 
filename_prov = '/Users/filipbirkfeldt/Desktop/StatStok/all_files/covProc.mat'; 
covProc = importdata(filename_prov); 
k_values=[1,2,3]; 
figure; 

for k=1:length(k_values)
    subplot(3,1,k_values(k)); 
    plot(covProc(1:end-k_values(k)), covProc(k_values(k)+1:end), '.');
    title(k_values(k)); 
end

% Estimation of the covariance funtion with xcov 
[ycov,lags]=xcov(covProc,20,'biased');
[ycov_norm, lags_norm] = xcov(covProc, 20, 'coeff'); 
% Cross-Variance : Hur covariancen varier mellan olika par av tiden. 
disp(ycov_norm)
% lengt(covProc) = 2000 , length(ycov_norm) = 41
% Hur de olika tidsparen varierar med varandra - om man följer t.ex t, t+1
% så syns den negativa covariansen tydligt. Följer man punkter i ycov_norm
% så syns det att covariansen är positiv.  

% *xcov calculates the covariance function for different values k, the 
% number of k:s tested is specified in the function. The covariance 
% function seems to have higher correlation between smaller 
% shifts in time rather than in larger. *


%% 3.3 Spectrum estimate of sum of harmonics 
% spekgui - estimate the covariance & spectral-densities 
% f_k = {5,10} sigma_k = {2,2}, phi e Rect(0,2pi), A e Rayleig(sigma^2)

% Q. Do the peaks have equal heigts? 
% Not equal heights, 

% function [rayamp]=enkelsumma(f,sigma2_,N,t,plotid)

f = [6 10];
sigma2_ = [2 2];
N = 1000;
dt = 1/(2*max(f)+1);
t = 0:dt:20;
enkelsumma(f,sigma2_,N,t,1);
%data.x = ans
%data.x = ans
spekgui
%data.x = rayamp(:,5);
%data.dt = dt;
% Q - väldigt lika i höjd

%Q. Draw rough sketches of the spectral density estimates obtained using 
% the periodogram. How do you explain the differences?

% En ny realisation ges vid varje anrop. För varje realisation så blir
% Olika frekvenser mer dominanta än de andra. 

% Q. The freq is not changing, how ever the phase and amplitude is drawn 
%from a set of possible outcomes. Which implies the spectral density of 
%the two functions may differ from each other. *

%% 4

filename_cello = '/Users/filipbirkfeldt/Desktop/StatStok/all_files/cello.mat'; 
cello = importdata(filename_cello); 
filename_tromb = '/Users/filipbirkfeldt/Desktop/StatStok/all_files/trombone.mat'; 
tromb = importdata(filename_tromb); 
 
% Q TROMBONEN
% 225, 445, 675 
% Antal övertoner på trombonen, 9st 

% Q. Can you see a strong noise peak at a particular frequency? 
% Ja, vid f=230

% CELLO

% Q. 
% f = 225, 445, 675, 900...
% Q 
% Ja, en fk = k*f0 
% Q 
% Om man kollar på den linjära skalan, tydlig topp vid f =900
% ---------------------------------
% *cello: 230, *
% - *trombone: 235, *
% 18 tones, 1 + 17, cello*
% - *6 tones, 1 + 5, Trombone*

% strong noise peeak - battery : Yes at 1/1.1? (trombone)* OR cello at 890 Hz ?, not really sure. *



%% 4.2
n = 2; 
cello2.x=cello.x(1:n:end);
cello2.dt=cello.dt*n;

%- Q. Has the spectrum changed? How has the spectrum range changed? At
% what frequencies do aliased peaks (if any) appear?
% *The range has ben reduced by half. *
% - *A peak has emerged between 1600 and 1700, and one between 1400 and 1500


% 1/7, n=7 ish 

% För trombonen 
n = 5; 
tromb2.x=tromb.x(1:n:end); 
tromb2.dt = tromb.dt*n; 

%Q, f=225, 450, 525, 675
%Q, när man samplar från n=9 och uppåt så får man inte
% med övertonerne längre. 
%%
% Down-sampla på rätt sätt istället: 
n=10;
soundsc(cello.x); 
cello2.x=decimate(cello.x,n);
cello2.dt=cello.dt*n; 
% Still aliased peaks (peak) bara väldigt låga frekvenser. endast 1 peak 
%%
soundsc(cello2.x) 

% Q. How much slower must this signal be sampled to give an aliasing in the spectrum?
% *The highest frequencies in the original cello sound was 4000 Hz, To 
% avoid aliasing you need to sample with a frequency of at least 2*4000, 
% so any sampling speed of below 8000 samples per second is too slow. *

% aliased peaks försvann 









        