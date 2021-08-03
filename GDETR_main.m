%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                      Replication Script for                     %%%%%
%%%%%                                                                 %%%%%
%%%%%     Growth, Distribution and Dynamic Inefficiency in Turkey     %%%%%
%%%%%     An Analysis of the Naïve Neoclassical Theory of Capital     %%%%%
%%%%%                                                                 %%%%%
%%%%%                        By  M. Aykut Attar                       %%%%%
%%%%%                                                                 %%%%%
%%%%%                          November 2020                          %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean:

clc; clear all; close all;

%% Load raw data:

D    = xlsread('GDETR_matlab.xls');

% The original data file organized and used by Altug et al. (2008, EREH) is
% located at 
%
%   http://myweb.sabanciuniv.edu/alpayf/files/2010/04/erehdata.xls
%
% and has 13 data columns. The GDETR_matlab.xls file this script calls at
% Line 7 should be organized in the following way:
% 
% (*) From Altug et al. (2008, EREH):
%
%       A1:A83 year              (A3:A85 in the original erehdata.xls file)
%       B1:B83 Real GDP          (B3:B85 in the original erehdata.xls file)
%       C1:C83 Capital stock     (E3:E85 in the original erehdata.xls file)
%       D1:D83 Employment        (H3:H85 in the original erehdata.xls file)
%       E1:E83 Human capital     (K3:K85 in the original erehdata.xls file)
% 
% (*) From the Penn World Tables v.9:
%
%       F1:F83 Human capital (PWT 9)
%              F1:F27  1923-1949  no data ( = NaN )
%              F28:F83 1950-2005  Turkey's "hc" data from the PWT 9

t    = D(:,1); % year
Y    = D(:,2); % Real GDP
K    = D(:,3); % Capital stock
L    = D(:,4); % Employment
h    = D(:,5); % Human capital (Altug et al. 2008)
hpwt = D(:,6); % Human capital (PWT 9)

%% Define output and capital per worker:

y    = Y./L;
k    = K./L;

%% Define the capital-output ratio:

kapa = K./Y;
[kapat,kapac] = hpfilter(kapa,100); % Hodrick-Prescott filter

%% Draw data figures:

figure(1)
subplot(2,1,1)
plot(t,log(Y),'or-')
hold on
plot(t,log(K),'sb-')
plot(t,log(L),'xk-')
hold off
xlim([1923 2005])
grid on
legend('Real GDP','Physical capital','Employment','location','NorthWest','orientation','vertical')
ylabel('logged levels')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];
subplot(2,1,2)
plot(t,log(y),'or-')
hold on
plot(t,log(k),'sb-')
plot(t,log(h),'xk-')
hold off
xlim([1923 2005])
ylim([-2 8.5])
grid on
legend('Real GDP per worker','Physical capital per worker','Human capital per worker','location','NorthWest','orientation','vertical')
ylabel('logged levels')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

figure(2)
plot(t,kapa,'or-')
hold on
plot(t,kapat,'k-')
plot([1923, 2005], [mean(kapa), mean(kapa)],'--','LineWidth',1,'Color','blue')
hold off
xlim([1923 2005])
grid on
legend('Capital-Output Ratio','H-P trend','sample average','location','NorthEast','orientation','vertical')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

figure(3)
plot(t,log(h),'or-')
hold on
plot(t,log(hpwt),'sb-')
hold off
xlim([1923 2005])
grid on
ylabel('logged levels')
legend('Human capital (Altug et al., 2008)','Human capital (PWT 9.0)','location','NorthWest','orientation','vertical')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

clear ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                              Growth                             %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate growth rates of y and k:

for z=1923:2004
   Gk(z-1922) = k(z-1922+1)/k(z-1922);
   Gy(z-1922) = y(z-1922+1)/y(z-1922);
end
Gk = Gk';
Gy = Gy';

%% Calculate average growth rates for 1923-29, 1930-49, 1950-79, 1980-2004:

Gyk = [Gy Gk];
aGyk = [mean(Gyk(1:7,:));mean(Gyk(8:27,:));mean(Gyk(28:57,:));mean(Gyk(58:82,:))];

%% Set alpha values: ( Y = K^[alpha] * (A h L)^[1-alpha] )

alfa1 = 1/3;  % Capital share     
alfa2 = 0.35;
alfa3 = 0.50;

%% Compute the TFP term A:

for z=1923:2005
   A1(z-1922) = ((Y(z-1922)/(K(z-1922)^alfa1))^(1/(1-alfa1)))*(1/(L(z-1922)));
   A1h(z-1922) = ((Y(z-1922)/(K(z-1922)^alfa1))^(1/(1-alfa1)))*(1/(L(z-1922)*h(z-1922)));
   A2(z-1922) = ((Y(z-1922)/(K(z-1922)^alfa2))^(1/(1-alfa2)))*(1/(L(z-1922)));
   A3(z-1922) = ((Y(z-1922)/(K(z-1922)^alfa3))^(1/(1-alfa3)))*(1/(L(z-1922)));
end
A1  = A1';   % Ignoring h, alfa = 1/3
A1h = A1h';
A2  = A2';   % Ignoring h, alfa = 0.35
A3  = A3';   % Ignoring h, alfa = 0.5

%% Compute the TFP term A with PWT human capital: 

for z=1923:2005
    if hpwt(z-1922) == NaN %#ok<*FNAN>
        A1pwt(z-1922) = NaN;
    else
        A1pwt(z-1922) = ((Y(z-1922)/(K(z-1922)^alfa1))^(1/(1-alfa1)))*(1/(L(z-1922)*hpwt(z-1922)));
    end
end
A1pwt  = A1pwt';

%% Draw TFP figures:

figure(4)
plot(t,log(A1),'or-')
hold on
plot(t,log(A1h),'sb-')
plot(t,log(A1pwt),'xk-')
hold off
xlim([1923 2005])
grid on
ylabel('logged levels')
legend('TFP','TFP corrected for h (Altug et al., 2008)','TFP corrected for h (PWT 9.0)','location','NorthWest','orientation','vertical')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

%% Growth accounting - Part 1: 1923-2005 A1, A2, A3

% Growth rates
for z=1923:2004
    GA1(z-1922) = A1(z-1922+1)/A1(z-1922);
    GA2(z-1922) = A2(z-1922+1)/A2(z-1922);
    GA3(z-1922) = A3(z-1922+1)/A3(z-1922);
end
GA1 = GA1';
GA2 = GA2';
GA3 = GA3';

% Average growth rates
aGA1 = [mean(GA1(1:7,:));mean(GA1(8:27,:));mean(GA1(28:57,:));mean(GA1(58:82,:))];
aGA2 = [mean(GA2(1:7,:));mean(GA2(8:27,:));mean(GA2(28:57,:));mean(GA2(58:82,:))];
aGA3 = [mean(GA3(1:7,:));mean(GA3(8:27,:));mean(GA3(28:57,:));mean(GA3(58:82,:))];

% Contributions to growth (%) 
Cont_A1 = 100*(((1-alfa1)*log(aGA1))./(log(aGyk(:,1))));
Cont_A2 = 100*(((1-alfa2)*log(aGA2))./(log(aGyk(:,1))));
Cont_A3 = 100*(((1-alfa3)*log(aGA3))./(log(aGyk(:,1))));

Cont_k1 = 100 - Cont_A1;
Cont_k2 = 100 - Cont_A2;
Cont_k3 = 100 - Cont_A3;

%% Growth accounting - Part 2: 1923-2005 A1h and h

% Growth rates
for z=1923:2004
    GA1h(z-1922) = A1h(z-1922+1)/A1h(z-1922);
    Gh(z-1922)   = h(z-1922+1)/h(z-1922);
end
GA1h = GA1h';
Gh   = Gh';

% Average growth rates
aGA1h = [mean(GA1h(1:7,:));mean(GA1h(8:27,:));mean(GA1h(28:57,:));mean(GA1h(58:82,:))];
aGh   = [mean(Gh(1:7,:));mean(Gh(8:27,:));mean(Gh(28:57,:));mean(Gh(58:82,:))];

% Contributions to growth (%)
Cont_A1h = 100*(((1-alfa1)*log(aGA1h))./(log(aGyk(:,1))));
Cont_h   = 100*(((1-alfa1)*log(aGh))./(log(aGyk(:,1))));

Cont_k1h = 100 - Cont_A1h - Cont_h;

%% Growth accounting - Part 3: 1950-2005 A1pwt hpwt

% Growth rates
for z=1923:2004
    if hpwt(z-1922) == NaN
        GA1pwt(z-1922) = NaN;
        Ghpwt(z-1922) = NaN;
    else
        GA1pwt(z-1922) = A1pwt(z-1922+1)/A1pwt(z-1922);
        Ghpwt(z-1922) = hpwt(z-1922+1)/hpwt(z-1922);
    end
end
GA1pwt = GA1pwt';
Ghpwt  = Ghpwt';

% Average growth rates
aGA1pwt = [mean(GA1pwt(1:7,:));mean(GA1pwt(8:27,:));mean(GA1pwt(28:57,:));mean(GA1pwt(58:82,:))];
aGhpwt  = [mean(Ghpwt(1:7,:));mean(Ghpwt(8:27,:));mean(Ghpwt(28:57,:));mean(Ghpwt(58:82,:))];

% Contributons to growth (%)
Cont_A1pwt = 100*(((1-alfa1)*log(aGA1pwt))./(log(aGyk(:,1))));
Cont_hpwt  = 100*(((1-alfa1)*log(aGhpwt))./(log(aGyk(:,1))));

Cont_k1hpwt = 100 - Cont_A1pwt - Cont_hpwt;

%% Display the growth accounting results:

clc;

disp('========================================================================');
disp('Table 1 - Growth Accounting (without corrections for human capital)     ');
disp('========================================================================');
disp('              Avg. Growth Rates (% per annum)    Contributions (%)      ');
disp('------------------------------------------------------------------------');
disp('Period        y       k       A                  k       A              ');
disp('------------------------------------------------------------------------');
disp('1923 - 1929   11.22   1.81    16.68              3.30    96.70          ');
disp('1930 - 1949    1.63  -0.03     2.98            -21.24   121.24          ');
disp('1950 - 1979    3.05   4.40     2.44             46.44    53.56          ');
disp('1980 - 2005    2.90   2.97     2.92             32.90    67.10          ');
disp('========================================================================');
disp('                                                                        ');
disp('                                                                        ');

disp('========================================================================');
disp('Table 2 - Growth Accounting (corrections for human capital)             ');
disp('========================================================================');
disp('              Avg. Growth Rates (% per annum)    Contributions (%)      ');
disp('------------------------------------------------------------------------');
disp('Period        y       k       A       h          k       A       h      ');
disp('------------------------------------------------------------------------');
disp('1923 - 1929   11.22   1.81    0.32    15.41      8.16    2.00    89.84  ');
disp('1930 - 1949    1.63  -0.03   -0.67     3.94    -31.64  -27.97   159.61  ');
disp('1950 - 1979    3.05   4.40   -0.33     2.79     46.03   -7.23    61.20  ');
disp('1980 - 2005    2.90   2.97    0.51     2.49     30.81   11.80    57.39  ');
disp('========================================================================');
disp('                                                                        ');
disp('                                                                        ');

disp('========================================================================');
disp('Table 3 - Growth Accounting (corrections for human capital, PWT 9.0)    ');
disp('========================================================================');
disp('              Avg. Growth Rates (% per annum)    Contributions (%)      ');
disp('------------------------------------------------------------------------');
disp('Period        y       k       A       h          k       A       h      ');
disp('------------------------------------------------------------------------');
disp('1950 - 1979   3.05    4.40    1.59    0.84       46.30   35.02   18.68  ');
disp('1980 - 2005   2.90    2.97    1.45    1.45       32.81   33.60   33.59  ');
disp('========================================================================');
disp('                                                                        ');
disp('                                                                        ');

%% Estimate long-run growth rates: ( ln(x_t) = ln(x_0) + t * ln[1 + g] + e_t )

% Altug et al. 2008 (1923-2005)
DepVar = [A1 A1h A2 A3 h y k Y];
tt = t-1922;
for i=1:8
    dv           = log(DepVar(:,i));
    Xv           = [ones(size(dv,1),1) tt(:,1)];
    b            = regress(dv,Xv);
    Gr(:,i)      = exp(b(2,1))-1;
    dvh          = b(1,1) + b(2,1)*[tt(:,1)];
    DepVarh(:,i) = dvh;
    clear dv Xv b dvh;
end
clear tt

% PWT (1950-2005)
DepVar_pwt = [A1pwt(28:83,1) hpwt(28:83,1)];
tt = 1:56; tt = tt';
for i=1:2
    dv               = log(DepVar_pwt(:,i));
    Xv               = [ones(size(dv,1),1) tt(:,1)];
    b                = regress(dv,Xv);
    Gr_pwt(:,i)      = exp(b(2,1))-1;
    dvh              = b(1,1) + b(2,1)*[tt(:,1)];
    DepVarh_pwt(:,i) = dvh;
    clear dv Xv b dvh;
end
clear tt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                           Distribution                          %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct factor prices:

for z = 1923:2005
   w1(z-1922) = (1-alfa1)*Y(z-1922)*(1/L(z-1922));
   r1(z-1922) = alfa1*Y(z-1922)*(1/K(z-1922));
end
w1 = w1';
r1 = r1';

for z = 1923:2005
   w2(z-1922) = (1-alfa2)*Y(z-1922)*(1/L(z-1922));
   r2(z-1922) = alfa2*Y(z-1922)*(1/K(z-1922));
end
w2 = w2';
r2 = r2';

for z = 1923:2005
   w3(z-1922) = (1-alfa3)*Y(z-1922)*(1/L(z-1922));
   r3(z-1922) = alfa3*Y(z-1922)*(1/K(z-1922));
end
w3 = w3';
r3 = r3';

%% Detrend the real wage: ( ws = w / Ah )

ws1 = w1./A1;
ws2 = w2./A2;
ws3 = w3./A3;

%% Construct the wage-rental ratio: ( wtor = ws / r )

wtor1 = ws1./r1;
[wtor1t,wtor1c] = hpfilter(wtor1,100);

wtor2 = ws2./r2;
[wtor2t,wtor2c] = hpfilter(wtor2,100);

wtor3 = ws3./r3;
[wtor3t,wtor3c] = hpfilter(wtor3,100);

%% Construct the Piketty differentials ( r - delt - gY )

for z = 1923:2004
    gY(z-1922) = (Y(z-1922+1)/Y(z-1922)) - 1;
end
gY = gY';

delt1 = 0.047;
delt2 = 0.042;
delt3 = 0.050;

PD1 = 100*(r1(1:82,1) - delt1 - gY);
[PD1t,PD1c] = hpfilter(PD1,100);
PD2 = 100*(r1(1:82,1) - delt2 - gY);
[PD2t,PD2c] = hpfilter(PD2,100);
PD3 = 100*(r1(1:82,1) - delt3 - gY);
[PD3t,PD3c] = hpfilter(PD3,100);

%% Draw the distribution figures

figure(5)
plot(t,wtor1,'or-')
hold on
plot(t,wtor1t,'k-')
plot([1923, 2005], [mean(wtor1), mean(wtor1)],'--','LineWidth',1,'Color','blue')
hold off
xlim([1923 2005])
grid on
legend('Wage-Rental Ratio','H-P trend','sample average','location','NorthWest','orientation','vertical')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

figure(6)
plot(t(1:82,1),PD1,'or-')
hold on
plot(t(1:82,1),PD1t,'k-')
plot([1923, 2004], [mean(PD1), mean(PD1)],'--','LineWidth',1,'Color','blue')
plot([1923, 2004], [0, 0],'LineWidth',2,'Color','green')
hold off
xlim([1923 2005])
grid on
ylabel('percent')
legend('PD (\delta = 4.7%)','H-P trend','sample average','location','SouthEast','orientation','vertical')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                        Dynamic Inefficiency                     %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct investment: ( I_t = K_t+1 - K_t + [delta] K_t )

for z = 1923:2004
   I1(z-1922) = K(z-1922+1) - K(z-1922) + delt1*K(z-1922);
   I2(z-1922) = K(z-1922+1) - K(z-1922) + delt2*K(z-1922);
   I3(z-1922) = K(z-1922+1) - K(z-1922) + delt3*K(z-1922);
end
I1  = I1';
I2  = I2';
I3  = I3';

%% Construct investment rate: ( s = I / Y )

for z = 1923:2004
   s1(z-1922) = I1(z-1922)/Y(z-1922);
   s2(z-1922) = I2(z-1922)/Y(z-1922);
   s3(z-1922) = I3(z-1922)/Y(z-1922);
end
s1  = s1';
s2  = s2';
s3  = s3';

[s1t,s1c] = hpfilter(s1,100);
[s2t,s2c] = hpfilter(s2,100);
[s3t,s3c] = hpfilter(s3,100);

%% Define capital income as profit: ( Pi = r * K )

Pi = r1.*K;

%% Construct Investment to Profit ratio:

I1toP = I1./Pi(1:82,:);
I2toP = I2./Pi(1:82,:);
I3toP = I3./Pi(1:82,:);

[I1toPt,I1toPc] = hpfilter(I1toP,100);
[I2toPt,I2toPc] = hpfilter(I2toP,100);
[I3toPt,I3toPc] = hpfilter(I3toP,100);

%% Draw the inefficiency figures

figure(7)
subplot(2,1,1)
plot(t(1:82,1),s1,'or-')
hold on
plot(t(1:82,1),s1t,'k-')
plot([1923, 2004], [mean(s1), mean(s1)],'--','LineWidth',1,'Color','blue')
plot([1923, 2004], [alfa1, alfa1],'LineWidth',2,'Color','green')
hold off
xlim([1923 2005])
grid on
legend('Inv-to-GDP (\delta = 4.7%)','H-P trend','sample average','location','NorthWest','orientation','horizontal')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];
subplot(2,1,2)
plot(t(1:82,1),I1toP,'or-')
hold on
plot(t(1:82,1),I1toPt,'k-')
plot([1923, 2004], [mean(I1toP), mean(I1toP)],'--','LineWidth',1,'Color','blue')
plot([1923, 2004], [1, 1],'LineWidth',2,'Color','green')
hold off
xlim([1923 2005])
grid on
legend('Inv-to-Profit (\delta = 4.7%)','H-P trend','sample average','location','NorthWest','orientation','horizontal')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                        Sensitivity Analyses                     %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8)
plot(t,log(A1),'or-')
hold on
plot(t,log(A2),'sb-')
plot(t,log(A3),'xk-')
hold off
xlim([1923 2005])
grid on
ylabel('logged levels')
legend('TFP (\alpha = 1/3)','TFP (\alpha = 0.35)','TFP (\alpha = 0.5)','location','SouthEast','orientation','vertical')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

figure(9)
subplot(2,1,1)
plot(t,wtor1,'or-')
hold on
plot(t,wtor2,'sb-')
plot(t,wtor3,'xk-')
hold off
xlim([1923 2005])
grid on
legend('Wage-Rental Ratio (\alpha = 1/3)','Wage-Rental Ratio (\alpha = 0.35)','Wage-Rental Ratio (\alpha = 0.5)','location','NorthWest','orientation','vertical')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];
subplot(2,1,2)
plot(t(1:82,1),PD1,'or-')
hold on
plot(t(1:82,1),PD2,'sb-')
plot(t(1:82,1),PD3,'xk-')
plot([1923, 2004], [0, 0],'LineWidth',2,'Color','green')
hold off
xlim([1923 2005])
grid on
ylabel('percent')
legend('PD (\delta = 4.7%)','PD (\delta = 4.2%)','PD (\delta = 5.0%)','location','SouthEast','orientation','vertical')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

figure(10)
subplot(2,1,1)
plot(t(1:82,1),s1,'or-')
hold on
plot(t(1:82,1),s2,'sb-')
plot(t(1:82,1),s3,'xk-')
plot([1923, 2004], [alfa1, alfa1],'LineWidth',2,'Color','green')
hold off
xlim([1923 2005])
grid on
legend('Inv-to-GDP (\delta = 4.7%)','Inv-to-GDP (\delta = 4.2%)','Inv-to-GDP (\delta = 5.0%)','location','NorthWest','orientation','horizontal')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];
subplot(2,1,2)
plot(t(1:82,1),I1toP,'or-')
hold on
plot(t(1:82,1),I2toP,'sb-')
plot(t(1:82,1),I3toP,'xk-')
plot([1923, 2004], [1, 1],'LineWidth',2,'Color','green')
hold off
xlim([1923 2005])
grid on
legend('Inv-to-Profit (\delta = 4.7%)','Inv-to-Profit (\delta = 4.2%)','Inv-to-Profit (\delta = 5.0%)','location','NorthWest','orientation','horizontal')
ax = gca;
ax.XTick = [1923 1930 1940 1950 1960 1970 1980 1990 2000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%