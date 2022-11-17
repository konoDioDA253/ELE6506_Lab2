% A =[3+4*j 0; 40*exp(j*deg2rad(50)) 6+5*j];
% [angle_A, rho_A] = cart2pol(real(A^2), imag(A^2));
% angle=rad2deg(angle_A)
% rho_A
% 

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
nom_image_coefs = ["2b_coefs_q0", "2b_coefs_q1", "2b_coefs_q2", "2b_coefs_q3", "2b_coefs_q4", "2b_coefs_q5"];
nom_image_ray=["2c_patron_q0", "2c_patron_q1", "2c_patron_q2", "2c_patron_q3", "2c_patron_q4", "2c_patron_q5"];
nom_image_dir="2d_directivite";
% Creation du repertoire des figures
mkdir('fig');

q=5;
D_max = zeros(q+1,1);

for q=0:1:q

    H=1e-2;
    PointNumber = 79;
    E0 = 1;
    lambda = 0.3;
    eta = 120*pi;
    beta = 2*pi/lambda;
    pas_theta = (pi/2)/(PointNumber-1);
    theta = 3*pi/4;
    X = zeros(PointNumber,1);
    Y = zeros(PointNumber,1);
    gamma = 0.5772;

    Z=zeros(PointNumber,PointNumber);
    I=zeros(PointNumber,1);
    V=zeros(PointNumber,1);
    E=zeros(PointNumber,1);
    Delta_c = lambda/10;
    Phi=zeros(PointNumber,1);

    pas_phiOrigine=360;
    Phi_origine = linspace(-pi/2,3*pi/2,pas_phiOrigine);
    E_phi = zeros(pas_phiOrigine,1);
    S=zeros(pas_phiOrigine,1);
    Smoy=0;
    D=zeros(pas_phiOrigine,1);


    % calculs de V et Z
    for i=1:1:PointNumber
        Phi(i)=theta;
        theta=theta+pas_theta;
        X(i) = 5*lambda*cos(theta) + 5*lambda/2;

        if i == 39
            Y(i) =0.030205579898166;
        else
            Y(i) = 5*lambda*sin(theta);
        end
    end

    %resolution du systeme
    for m=1:1:PointNumber
        Rho = sqrt(X(m)^2+Y(m)^2);
        if (Y(m)==0)
            Phi(m)=pi;
        end
        if (Y(m)<0)
            Phi(m)=pi+atan(abs(Y(m)/X(m)));
        end
        if (Y(m)>0)
            Phi(m)=pi-atan(abs(Y(m)/X(m)));
        end
        V(m)=E0*(cos(Phi(m))^q)*exp(-j*beta*Rho)/sqrt(Rho);        



        for n=1:1:PointNumber        
            if (m==n)
                Z(m,n) = (beta*eta*Delta_c/4) * (1-j*(2/pi) * log(gamma*beta*Delta_c/(4*exp(1))) );
            else
                argument = beta*sqrt((X(n)-X(m))^2 + (Y(n)-Y(m))^2);
                Z(m,n) = (beta*eta*Delta_c/4) * besselh(0,2,argument) ;
            end
        end
    end

    I=Z\V;


    %Calcul de E et S
    for i=1:1:pas_phiOrigine
        if ((Phi_origine(i) <= pi/2) &&( Phi_origine(i)>=(-pi/2)))
            Esource = 0;
        end
        if ((Phi_origine(i) >= pi/2) && (Phi_origine(i)<=3*pi/2))
%             R=1000*cos(Phi_origine(i))+1000*sin(Phi_origine(i));
            Esource = E0*(cos(Phi_origine(i)))^q;
        end
    %     S_source=0.5*real((Esource)/sqrt(R))^2/eta);    
        E_Total_tubulaire=0;
        for n=1:1:PointNumber
            R=X(n)*cos(Phi_origine(i))+Y(n)*sin(Phi_origine(i));
            E_Total_tubulaire=E_Total_tubulaire-eta*I(n)*exp(j*beta*R)*sqrt(j*beta/(8*pi));%/sqrt(R);
    %         S(n)=0.5*real((-eta*I(n)*exp(-j*beta*R)*sqrt(j*beta/(8*pi))/sqrt(R))^2/eta);
    %         if n==38
    %             disp("here!")
    %         end
    %         E_Total_tubulaire=E_Total_tubulaire+j*eta*I(n)*1*(cos(beta*H*cos(Phi_origine(i))/2)-cos(beta*H/2))/(2*pi*(H/2)*sin(beta)*sin(Phi_origine(i)));
        end
        E_phi(i) = Esource+E_Total_tubulaire;
        S(i)=0.5*real(E_phi(i).^2)/eta;
        Smoy=Smoy+H*S(i)*pas_phiOrigine;
    end
    Smoy=Smoy/(2*pi);

    E_phi_norme = abs(E_phi);
%     D=S./Smoy;
    Emoy_norme = abs(mean(E_phi,'all'));
    D=E_phi_norme.^2./Emoy_norme.^2;
    
    D_max(q+1) = max(D);
    Module_I = abs(I);
    Angle_I = rad2deg(angle(I));
    Phi = rad2deg(Phi);

    %% 2 graphiques en subplot
    fig = figure();
    subplot(2,1,1);
    plot(Phi,Module_I, 'b-o') % Trace de la reference
    grid minor % Jolie grille
    string_titre = strcat('$|I_n|$ en fonction de $\phi$ pour q = ', num2str(q));
    string_Y = 'coefficients $|I_n|$';
    title(string_titre,'Interpreter','latex','FontSize',16) % On veut en LaTeX le titre et en 16pt
    xlabel('$\phi$  ($^\circ$)','Interpreter','latex','FontSize',12) % abscisses
    ylabel(string_Y,'Interpreter','latex','FontSize',12) % Ordonnee
    % Legend, plus delicat, il faut une structure et placer cette derniere (on
    % choisit l'option ?best? de placement automatique. En plus il y a les
    % accents!
    % string = "\Phi(t) bruit?"; % Deux entrees pour la legende (dans l'ordre!)
    legend(string_titre,'Interpreter','latex','FontSize',12,'Location','best')

    subplot(2,1,2); % Trace pour la commande, pas les memes unites ! Jamais sur le meme graphique!!
    grid minor
    hold on
    string_titre = strcat('$\angle I_n$ en fonction de $\phi$ pour q = ', num2str(q));
    string_Y = '$\angle I_n$ ($^\circ$)';
    plot(Phi,Angle_I, 'b-o') % Trace de la commande
    title(string_titre,'Interpreter','latex','FontSize',16)
    xlabel(' $\phi$ ($^\circ$)','Interpreter','latex','FontSize',12) 
    ylabel(string_Y,'Interpreter','latex','FontSize',12) 
    legend(string_titre,'Interpreter','latex','FontSize',12,'Location','best')

    fig2 = figure();
    polarplot(Phi_origine,E_phi_norme)
    string_titre = strcat('Patron de rayonnement en fonction de $\phi$ pour q = ', num2str(q));
    title(string_titre,'Interpreter','latex','FontSize',16) % On veut en LaTeX le titre et en 16pt

    % fig3 = figure();
    % grid minor
    % hold on
    % string_titre = strcat('Directivite en fonction de q pour q = ', num2str(q));
    % string_Y = ' ($^\circ$)';
    % plot(q,D_max, 'b-o') % Trace de la commande
    % title(string_titre,'Interpreter','latex','FontSize',16)
    % xlabel(' $\phi$ ($^\circ$)','Interpreter','latex','FontSize',12) 
    % ylabel(string_Y,'Interpreter','latex','FontSize',12) 
    % % Pas besoin de stucture ici, une seule legende
    % legend(string_titre,'Interpreter','latex','FontSize',12,'Location','best')


    %% Section 3 -- Export de la figure (automatique, eh oui!)
    %1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
    print(fig,'-depsc',strcat(path,nom_image_coefs(q+1),'.eps'))
    %2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
    print(fig,'-dpng','-r600',strcat(path,nom_image_coefs(q+1),'.png'))

    %1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
    print(fig2,'-depsc',strcat(path,nom_image_ray(q+1),'.eps'))
    %2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
    print(fig2,'-dpng','-r600',strcat(path,nom_image_ray(q+1),'.png'))

    % %1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
    % print(fig3,'-depsc',strcat(path,nom_image_dir(q+1),'.eps'))
    % %2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
    % print(fig3,'-dpng','-r600',strcat(path,nom_image_dir(q+1),'.png'))

end

q=0:1:5;
fig3 = figure();
grid minor
hold on
string_titre = 'Directivite en fonction de q';
string_Y = ' Directivite';
plot(q,D_max, 'b-o') % Trace de la commande
title(string_titre,'Interpreter','latex','FontSize',16)
xlabel(' Valeur de q (entier)','Interpreter','latex','FontSize',12) 
ylabel(string_Y,'Interpreter','latex','FontSize',12) 
% Pas besoin de stucture ici, une seule legende
legend(string_titre,'Interpreter','latex','FontSize',12,'Location','best')

%1) en format vectoriel .eps pour les vrais qui utilisent LaTeX
print(fig3,'-depsc',strcat(path,nom_image_dir,'.eps'))
%2) en format PNG (haute resolution 600dpi, pour ceux qui utilisent Word)
print(fig3,'-dpng','-r600',strcat(path,nom_image_dir,'.png'))






