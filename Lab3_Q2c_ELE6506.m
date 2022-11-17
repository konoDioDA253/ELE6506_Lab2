clear all
close all
clc

% Choix du systeme d'exploitation
is_unix = 0; % 0 si votre systeme est windows; 1 sinon (Mac,Ubuntu)
if is_unix
    % Unix
    path = "fig/";
else
    % Windows
    path = "fig\";
end

% Creation des noms des figures
nom_image1="1a";
nom_image2="1c";
% nom_image3="2a_F_Fd";
% Creation du repertoire des figures
mkdir('fig');

% psi = tf('s');
% syms psi;
% AF_norm=sin(N*psi/2)/(N*sin(psi/2));
% f = 20*log10(AF_norm)+11;
% zero=fzero(f, 0);
d_sur_lambda = 0.5;

psi_11 = zeros(6,3);
cos_theta_min = zeros(6,3);
theta_min = zeros(6,3);
len=1000;
theta = linspace(0,360,len);
global N;
for N=5:1:10
    for i=1:1:3
        d_sur_lambda = 0.4 + i*0.1;
        fun = @f; % function
        x0 = 0.5; % initial point
        psi_11(N-4,i) = fzero(fun,x0);
        cos_theta_min(N-4,i) = (2*pi-psi_11(N-4,i))/(2*pi*d_sur_lambda)-1;
        theta_min(N-4,i) = rad2deg(acos(cos_theta_min(N-4,i)));
        if (i == 3) && (N == 10)
            d_sur_lambda = 0.7;
            alpha = 2*pi-2*pi*d_sur_lambda*cos_theta_min(6,3);
            psi_courant = alpha+2*pi*d_sur_lambda*cos(deg2rad(theta));
            AF=abs(sin(N*psi_courant/2)./sin(psi_courant/2));
        end
    end
%     fun = @f; % function
%     x0 = 0.5; % initial point
%     psi_11(N-4) = fzero(fun,x0);
%     cos_theta_min(N-4) = (2*pi-psi_11(N-4))/(2*pi*d_sur_lambda)-1;
%     theta_min(N-4) = rad2deg(acos(cos_theta_min(N-4)));
end


N=5:1:10;

fig = figure();
% subplot(2,1,1);
plot(N,theta_min(:,1),'b-o') % Trace de la reference
hold on
plot(N,theta_min(:,2),'r-o') % Trace de la reference
plot(N,theta_min(:,3),'m-o') % Trace de la reference
grid minor % Jolie grille
string_titre ='$\theta_{min}$ en fonction de $N$ pour differentes valeurs de $d/\lambda$';
string_Y = '$\theta_{min}$';
title(string_titre,'Interpreter','latex','FontSize',16) % On veut en LaTeX le titre et en 16pt
xlabel('$N$','Interpreter','latex','FontSize',12) % abscisses
ylabel(string_Y,'Interpreter','latex','FontSize',12) % Ordonnee
% Legend, plus delicat, il faut une structure et placer cette derniere (on
% choisit l'option ?best? de placement automatique. En plus il y a les
% accents!
% string = "\Phi(t) bruit?"; % Deux entrees pour la legende (dans l'ordre!)
legend('d/$\lambda$ = 0.5','d/$\lambda$ = 0.6','d/$\lambda$ = 0.7','Interpreter','latex','FontSize',12,'Location','best')

fig2 = figure();
polarplot(theta,AF)
string_titre = 'Facteur de reseau pour N=10 et d/$\lambda$ = 0.7';
title(string_titre,'Interpreter','latex','FontSize',16) % On veut en LaTeX le titre et en 16pt

    %% Section 3 -- Export de la figure (automatique, eh oui!)
%1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
print(fig,'-depsc',strcat(path,nom_image1,'.eps'))
%2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
print(fig,'-dpng','-r600',strcat(path,nom_image1,'.png'))

%1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
print(fig2,'-depsc',strcat(path,nom_image2,'.eps'))
%2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
print(fig2,'-dpng','-r600',strcat(path,nom_image2,'.png'))














%% Question 2
clear all
close all

% Choix du systeme d'exploitation
is_unix = 0; % 0 si votre systeme est windows; 1 sinon (Mac,Ubuntu)
if is_unix
    % Unix
    path = "fig/";
else
    % Windows
    path = "fig\";
end

% Creation des noms des figures
nom_image1="2c_wn_an";
nom_image2="2c_im";
nom_image3="2c_F_Fd";
% Creation du repertoire des figures
mkdir('fig');


f = 1.9e9;
c=3e8;
lambda = c/f;
L=2;
d=0.6*lambda;
w_n=zeros(15,1);
a_n = zeros(15,1);
n=1;
max=0;
w_n(1)=1/(L/lambda);
P=25;
while(1)
    w_n(n+1)=(n+1)/(L/lambda);
%     if (( cos(135*pi/180) < w_n(n+1)) || (w_n(n+1) < cos(95*pi/180)) )
%         a_n(n) = 1/w_n(n);
%     else
%         a_n(n)=0;
%     end
%     i_m = ;
    n=n+1;
    if (n>max)
        max = n;
    end
    
    if ~(abs(w_n(n)) <= 1)
        break;
    end
end


max = max-1;
w_n = w_n(1:12);
a_n = a_n(1:12);
i_m = zeros(25,1);
z_m = zeros(25,1);
M=-max:1:max;
W_n = [-flipud(w_n); 0; w_n];
% A_n = [-flipud(a_n); 0; a_n];

%%%
for n=1:1:25
    if ( cos(135*pi/180) <= W_n(n)) && (W_n(n) <= cos(95*pi/180)) 
        denominateur = W_n(n)*sin(acos(deg2rad(W_n(n))));
        a_n(n) = abs(1/denominateur);
    else
        a_n(n)=0;
    end
%     n=n+1;    
%     if ~(abs(W_n(n)) <= 1)
%         break;
%     end
end


somme=0;
for i = 1:1:25
    m = M(i);
    z_m(i) = m*d;
    somme=0;
    for k = 1:1:25
        n=M(k);
        somme=somme + a_n(k)*exp(-1j*2*pi*(z_m(i)/lambda)*W_n(k));
    end
    i_m(i)=somme/P;
end

len=1000;
somme = 0;
theta = linspace(0,180,len);
w=zeros(len,1);
F=zeros(len,1);
Fd=zeros(len,1);
for index=1:1:len
    w(index)=cos(deg2rad(theta(index)));
    somme = 0;
    for k=1:1:25
        n=M(k);
        somme = somme + a_n(k)*sin(0.5*P*(w(index)-W_n(k)) * 2*(pi*d/lambda))/(P*12.07*sin(0.5*((w(index)-W_n(k))  * (2*pi*d/lambda) )));
    end
    F(index)=abs(somme);
    if (theta(index) < 135) && (theta(index) > 95)
        Fd(index)=abs(1/(12.07*cos(deg2rad(theta(index))*sin(deg2rad(theta(index)))) ));
    else
        Fd(index)=0;
    end
end

%% 2 graphiques en subplot
fig = figure();
subplot(2,1,1);
plot(-12:1:12,W_n, 'b-o') % Trace de la reference
grid minor % Jolie grille
string_titre ='$w_n$ en fonction de n';
string_Y = 'coefficients $w_n$';
title(string_titre,'Interpreter','latex','FontSize',16) % On veut en LaTeX le titre et en 16pt
xlabel('n','Interpreter','latex','FontSize',12) % abscisses
ylabel(string_Y,'Interpreter','latex','FontSize',12) % Ordonnee
% Legend, plus delicat, il faut une structure et placer cette derniere (on
% choisit l'option ?best? de placement automatique. En plus il y a les
% accents!
% string = "\Phi(t) bruit?"; % Deux entrees pour la legende (dans l'ordre!)
legend(string_titre,'Interpreter','latex','FontSize',12,'Location','best')

subplot(2,1,2); % Trace pour la commande, pas les memes unites ! Jamais sur le meme graphique!!
grid minor
hold on
string_titre = '$a_n$ en fonction de n';
string_Y = 'coefficients $a_n$';
plot(-12:1:12,a_n, 'b-o') % Trace de la commande
title(string_titre,'Interpreter','latex','FontSize',16)
xlabel('n','Interpreter','latex','FontSize',12) % abscisses
ylabel(string_Y,'Interpreter','latex','FontSize',12) 
legend(string_titre,'Interpreter','latex','FontSize',12,'Location','best')


fig2 = figure();
subplot(2,1,1);
plot(-12:1:12,abs(i_m), 'b-o') % Trace de la reference
grid minor % Jolie grille
string_titre ='$|i_m|$ en fonction de m';
string_Y = 'coefficients $|i_m|$';
title(string_titre,'Interpreter','latex','FontSize',16) % On veut en LaTeX le titre et en 16pt
xlabel('n','Interpreter','latex','FontSize',12) % abscisses
ylabel(string_Y,'Interpreter','latex','FontSize',12) % Ordonnee
% Legend, plus delicat, il faut une structure et placer cette derniere (on
% choisit l'option ?best? de placement automatique. En plus il y a les
% accents!
% string = "\Phi(t) bruit?"; % Deux entrees pour la legende (dans l'ordre!)
legend(string_titre,'Interpreter','latex','FontSize',12,'Location','best')

subplot(2,1,2); % Trace pour la commande, pas les memes unites ! Jamais sur le meme graphique!!
grid minor
hold on
string_titre = '$\angle i_m$ en fonction de m';
string_Y = '$\angle i_m$';
plot(-12:1:12,rad2deg(angle(i_m)), 'b-o') % Trace de la commande
title(string_titre,'Interpreter','latex','FontSize',16)
xlabel('m','Interpreter','latex','FontSize',12) % abscisses
ylabel(string_Y,'Interpreter','latex','FontSize',12) 
legend(string_titre,'Interpreter','latex','FontSize',12,'Location','best')


fig3 = figure();
% subplot(2,1,1);
plot(theta,F,'b-') % Trace de la reference
hold on
plot(theta,Fd,'r--') % Trace de la reference
grid minor % Jolie grille
string_titre ='F et Fd en fonction de $\theta$';
string_Y = 'F et Fd';
title(string_titre,'Interpreter','latex','FontSize',16) % On veut en LaTeX le titre et en 16pt
xlabel('$\theta$','Interpreter','latex','FontSize',12) % abscisses
ylabel(string_Y,'Interpreter','latex','FontSize',12) % Ordonnee
% Legend, plus delicat, il faut une structure et placer cette derniere (on
% choisit l'option ?best? de placement automatique. En plus il y a les
% accents!
% string = "\Phi(t) bruit?"; % Deux entrees pour la legende (dans l'ordre!)
legend('F','Fd','Interpreter','latex','FontSize',12,'Location','best')

%% Section 3 -- Export de la figure (automatique, eh oui!)
%1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
print(fig,'-depsc',strcat(path,nom_image1,'.eps'))
%2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
print(fig,'-dpng','-r600',strcat(path,nom_image1,'.png'))

%1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
print(fig2,'-depsc',strcat(path,nom_image2,'.eps'))
%2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
print(fig2,'-dpng','-r600',strcat(path,nom_image2,'.png'))

%1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
print(fig3,'-depsc',strcat(path,nom_image3,'.eps'))
%2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
print(fig3,'-dpng','-r600',strcat(path,nom_image3,'.png'))

