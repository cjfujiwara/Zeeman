function shift_register

% Lowest energies of -9,-7,-5
B=200.02;
D3= -760.5844;
D2= -807.4460;
D1= -851.9368;

% Subtract energy
D1r=D1-D1;
D2r=D2-D1;
D3r=D3-D1;


%% Static matrices
H0=2*pi*[D3r 0 0;
    0 D2r 0;
    0 0 D1r];

Ix=[0 1 0;
    1 0 1;
    0 1 0];


%% Define Perturbation

omega_lim=2*pi*[43 48];
omega_vec=linspace(omega_lim(1),omega_lim(2),1E3);

% Rabi frequency
Omega=2*pi*0.5;

Eq=zeros(length(omega_vec),3);

for kk=1:length(omega_vec)
    omega=omega_vec(kk);
    Eq(kk,1)=H0(1,1)-2*omega;
    Eq(kk,2)=H0(2,2)-1*omega;
end

%% Plot

figure(10);
clf
plot(omega_vec/(2*pi),Eq(:,1)/(2*pi),'linewidth',1)
hold on
plot(omega_vec/(2*pi),Eq(:,2)/(2*pi),'linewidth',1)
plot(omega_vec/(2*pi),Eq(:,3)/(2*pi),'linewidth',1)

legend({'E_{-5/2} - 2h\nu','E_{-7/2} - h\nu','E_{-9/2}'})

xlabel('frequency (MHz)');
ylabel('energy (MHz)');


%% Floquet
vFAll=zeros(3,3,length(omega_vec));

% number of drive steps
n=20;
thetaVec=linspace(0,2*pi,n);

UAll=ones(3,3,length(omega_vec));
    
for nn=1:length(omega_vec) 
    % Time perturbing hamiltonian
    H1=@(theta) Ix*Omega*cos(theta);

    omega=omega_vec(nn);      
    Tau=2*pi/omega;
    dTau=Tau/n;

    for kk=1:length(thetaVec)
       Htot=H0+H1(thetaVec(kk));
       Unow=expm(-1i*Htot*dTau);       
       UAll(:,:,kk)=Unow;    
    end

    UCycle=eye(3);
    for kk=1:n
        UCycle=UAll(:,:,kk)*UCycle;
    end

    G=1i*logm(UCycle)/Tau;
    [vF,b]=eig(full(G));
    Efloquet=real(b)*ones(size(vF,1),1);

    
%     % Sort by probability overlap with original states
    vS=eye(3);              % static eigenvectors
    % Compare floquet to statis spectrum
    Cvec=vS*conj(vF);                       % Compute overlap
    Pvec=Cvec.*conj(Cvec);                  % Find probability          
    %Acquire the highest probability floquet states
    [Pmax,inds]=max(Pvec,[],2);             % Sort by overlap with static
%     
% 
%     % sort by symmetry vs antisymmetry
%     if nn>1
%         vF(:,1)=vF(:,1)*sign(real(vF(2,1)));
%         vF(:,2)=vF(:,2)*sign(real(vF(2,2)));   
%     end
%     S=sign(real(vF));
%     P=prod(S,1);
%     [~,inds]=sort(P,2);
    
 
    
%     Efloquet=Efloquet(inds);               % Sort and  energies      
%     vF=vF(:,inds);    
    vFAll(:,:,nn)=vF;    
    engs(:,nn)=Efloquet;
    
end

strH='$H=2\pi\pmatrix{f_{59} & 0 & 0 \cr 0 & f_{79} & 0 \cr 0 &0 &0}+\Omega\cos(2\pi f t)\pmatrix{0 & 1 & 0 \cr 1 & 0 & 1 \cr 0 & 1 & 0}$';
strP=['$\Omega=2\pi\times' num2str(Omega/(2*pi)) '~\mathrm{MHz},' ...
    '~B=200~\mathrm{G}$'];


hf=figure(11);

hf.Color='w';
clf

co=get(gca,'colororder');
title('rf floquet analysis');

plot(omega_vec/(2*pi),engs(1,:)/(2*pi),'k.');
hold on
plot(omega_vec/(2*pi),engs(2,:)/(2*pi),'k.');
p4=plot(omega_vec/(2*pi),engs(3,:)/(2*pi),'k.');

p1=plot(omega_vec/(2*pi),Eq(:,1)/(2*pi),'linewidth',2,'color',co(1,:));
hold on
p2=plot(omega_vec/(2*pi),Eq(:,2)/(2*pi),'linewidth',2,'color',co(2,:));
p3=plot(omega_vec/(2*pi),Eq(:,3)/(2*pi),'linewidth',2,'color',co(3,:));


xlabel('frequency (MHz)');
ylabel('quasienergy (MHz)');

set(gca,'fontsize',14,'fontname','times');

legend([p1 p2 p3 p4],{'$\epsilon_{-5/2} - 2h\nu$','$\epsilon_{-7/2} - h\nu$','$\epsilon_{-9/2}$','$\widetilde{\epsilon}$'},...
    'interpreter','latex')

text(0.02,0.02,strH,'units','normalized','fontsize',12,'interpreter','latex',...
    'verticalalignment','bottom');
text(0.02,0.98,strP,'units','normalized','fontsize',12,'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left');

end

