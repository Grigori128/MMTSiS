%% Zadanie 1

clear variables
syms s12 s13 s23

%przypadek 2x2
A = [0 1; -1 -2];
Q = eye(2);
S = [0 s12; -s12 0];

L = (A')*S + S*A;
R = (A')*Q - Q*A;

s12 = 1;
S = [0 s12; -s12 0];
P = 0.5*(S-Q)*inv(A);
det(P);

%sprawdzenie
P1 = lyap(A',Q);

%% Zadanie 2 Obliczanie L

clear variables
syms y1 y2 real

%definiowanie macierzy
A = [0 1; 4 -1];
B = [0; 1];
Y = [y1 y2];

%wybór macierzy Q (dodatnio okreœlona i symetryczna)
Q = [1 -1; -1 2];

%obliczenie macierzy L
L = Q*(A')+ A*Q + (Y')*(B') + B*Y;

%% Zadanie 2 Wykres Y z obszarem rozwi¹zañ


y1 = -20:0.1:20;
y2 = -20:0.1:20;
[Y1,Y2] = meshgrid(y1,y2);
detL = -Y1.^2 - 14*Y1 - 4*Y2 - 25;

figure(1)
hold on
mesh(Y1,Y2,detL)
patch([20 -20 -20 20],[20 20 -20 -20],'r','FaceAlpha',0.2)
plot3(0,-7,3,'b+','MarkerSize',6,'LineWidth',2);
view(39,37);
xlabel('y_1')
ylabel('y_2')
zlabel('det(L)')
hold off
grid on;


%% Zadanie 2 Obliczanie K i Ap

%wybór punktu (y1,y2) z obszaru rozwi¹zañ, wyliczenie K
y1 = 0;
y2 = -7;
K = [y1 y2]*inv(Q);
k1 = K(1);
k2 = K(2);

%obliczenie nowej macierzy g³ównej uk³adu
Ap = A + [ 0 0; k1 k2];
%obliczenie wartoœci w³asnych
lambda=eig(Ap);


%% Zadanie 2 Wykresy z Simulinka
out = sim('Lapsim2');

figure(2)
subplot(2,1,1)
plot(out.ScopeData.time,out.ScopeData.signals(1,1).values);
ylim([-0.1 0.4])
xlabel('t (s)');
ylabel('x_1');
grid on;

subplot(2,1,2)
plot(out.ScopeData.time,out.ScopeData.signals(1,2).values);
ylim([-0.1 0.2])
xlabel('t (s)');
ylabel('x_2');
grid on;

%% Zadanie 3 Obliczanie L

clear variables
syms s12 s13 s23 y1 y2 y3 real

%definiowanie macierzy
A = [0 1 0 ; 0 -0.875 -20 ; 0 0 -50];
B = [0 ; 0 ; 50];
Y = [y1 y2 y3];

%sprawdzenie stabilnoœci - det(A)=0 wiêc wartoœci w³asne
lambda0 = eig(A);

%wybór macierzy Q (dodatnio okreœlona i symetryczna)
Q = [2 -3 0; -3 6 2; 0 2 3]; %kluczowe q23

%sprawdzenie dodatniej okreœlonoœci Q
display(det(Q(1,1)))
display(det(Q(1:2,1:2)))
display(det(Q))

%obliczenie macierzy L
L = Q*(A') + A*Q + (Y')*(B') + B*Y;

%sprawdzenie ujemnej okreœlonoœci L
display(det(L(1,1)))
display(det(L(1:2,1:2)))
display(det(L))

%nierównoœæ det(L) < 0
%226250*r1^2 + 43125*r1*r2 - (971275*r1)/8 + 15000*r2^2 - 95325*r2 +(749775*r3)/16 + 178835/16 < 0


%% Zadanie 3 Obliczanie K i Ap

%wybór punktu (y1,y2,y3) z obszaru rozwi¹zañ, wyliczenie K
y1 = 0;
y2 = 0;
y3 = -1;
K = [y1 y2 y3]*inv(Q);
k1 = K(1);
k2 = K(2);
k3 = K(3);

%obliczenie nowej macierzy g³ównej uk³adu
Ap = A + B*[k1 k2 k3];
%obliczenie wartoœci w³asnych
lambda=eig(Ap);


%% Zadanie 3 Sprawdzenie stabilnoœci uk³adu ze sprzê¿eniem

%zdefiniowanie macierzy skoœnie symetrycznej
S = [0 s12 s13; -s12 0 s23; -s13 -s23 0];

L = (Ap')*S + S*Ap;
R = (Ap')*Q - Q*Ap;

%obliczenie macierzy S
s = inv([-7/8 200 -300; -20 -200 0; 0 1 -1607/8])*[595.3750; 840; 1.11825e+03];
s12 = s(1);
s13 = s(2);
s23 = s(3);
S = [0 s12 s13; -s12 0 s23; -s13 -s23 0];

%obliczenie macierzy P
P = 0.5*(S-Q)*inv(Ap);
%sprawdzenie minorów
M1 = det(P(1,1));
M2 = det(P(1:2,1:2));
M3 = det(P);

%% Zadanie 3 Wykresy z Simulinka

t = 10;
out = sim('Lapsim3');

figure(4)
subplot(3,1,1)
plot(out.ScopeData.time,out.ScopeData.signals(1,1).values);
xlabel('t (s)');
ylabel('\theta (rad)');
set(gca,'LooseInset',get(gca,'TightInset'));
grid on;

subplot(3,1,2)
plot(out.ScopeData.time,out.ScopeData.signals(1,2).values);
xlabel('t (s)');
ylabel('\omega (rad/s)');
set(gca,'LooseInset',get(gca,'TightInset'));
grid on;

subplot(3,1,3)
plot(out.ScopeData.time,out.ScopeData.signals(1,3).values);
xlabel('t (s)'); 
ylabel('\tau (Nm)');
set(gca,'LooseInset',get(gca,'TightInset'));
grid on;

%print(1, '-dpng','odpowiedzi3', '-r800')


%% Zadanie 3 Model F15 w STL
fv = stlread('F-15_jet_plane.stl');
%fv.vertices = fv.vertices*[cos(pi),0,sin(pi);0 1 0; -sin(pi) 0, cos(pi)];
ftmp = fv;

%% Zadanie 3 Wizualizacja F15 i symulacja ruchu

figure(5)
set(5,'Position', [200 200 1300 600]);
set(gca, 'LooseInset', [0,0,0,0]);
grid on;
time = out.ScopeData.time;
theta = out.ScopeData.signals(1,1).values;
omega = out.ScopeData.signals(1,2).values;
tau = out.ScopeData.signals(1,3).values;
ctrl = out.ScopeData.signals(1,4).values;

for r = 1:5:length(omega)
set(0,'CurrentFigure',5)
title(strcat('t =',{' '},string((r*0.01)-0.01),'s'));

subplot(4,2,1)
plot(time(1:r),theta(1:r),'k-','LineWidth',1.8);
ylim([min(theta)-0.1*max(theta) max(theta)+0.1*max(theta)])
xlim([0 t]);
xlabel('t (s)');
ylabel('\theta (rad)');
grid on;

subplot(4,2,3)
plot(time(1:r),omega(1:r),'b-','LineWidth',1.8);
ylim([min(omega)-0.1*max(omega) max(omega)+0.1*max(omega)])
xlim([0 t]);
xlabel('t (s)');
ylabel('\omega (rad/s)');
grid on;

subplot(4,2,5)
plot(time(1:r),tau(1:r),'r-','LineWidth',1.8);
ylim([min(tau)-0.1*max(tau) max(tau)+0.1*max(tau)])
xlim([0 t]);
xlabel('t (s)'); 
ylabel('\tau (Nm)');
grid on;

subplot(4,2,7)
plot(time(1:r),ctrl(1:r),'g-','LineWidth',1.8);
xlim([0 t]);
ylim([min(ctrl)-(0.1*max(ctrl)+0.1) max(ctrl)+(0.1*max(ctrl)+0.1)])
xlabel('t (s)'); 
ylabel('r');
grid on;

subplot(4,2,[2,4,6,8])
Rz = [cos(theta(r)),0,sin(theta(r));0 1 0; -sin(theta(r)) 0, cos(theta(r))];
cla
grid on;
view(170,4); % <--- blokada kamery
ftmp.vertices = fv.vertices*Rz;
axis([-8 8 -3 20 -10 10]);
patch([20 -20 -20 20],[20 20 -20 -20],'k','FaceAlpha',0.2)
camlight('left');
material('metal');
patch(ftmp,'FaceColor',       [0 1 0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
material('metal');
drawnow
%drawnow limitrate
%pause(0.001)
% if r == 150
% print(1, '-dpng','wiz', '-r800')
% end
end

