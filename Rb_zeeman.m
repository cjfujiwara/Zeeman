function [out,hF1,hF2,hF3]=Rb_zeeman
% Author : C Fujiwara
%
% Calculate Zeeman shifts for 87Rb
%
% Output the eigen energies for a variety of magnetic fields Bvec. The
% output structure "out" is assigned the eigen values
out=struct;
Bvec=linspace(0,15000,1E5);
out.B=Bvec;


%% Constants
% Bohr magneton in MHz/Gauss
muB=1.39962449 ;    

% 40K has S=1/2 and I=4;
S=0.5;
I=3/2;

%% Coupling Constants
a_s12=3.41734130545215*1E3;     % 5S_1/2 magnetic dipole hyperfine
a_p12=408.328;                  % 5P_1/2 magnetic dipole hyperfine
a_p32=84.7185;                  % 5P_3/2 magnetic dipole hyperfine
b_p32=12.4965;                  % 5P_3/2 electric quadrupole hyperfine

gS=2.0023193043622; % Electron gyromagnetic ratio
gL=1;               % Electron orbital gyromagnetic ratio
gI=-0.0009951414;     % 40K nucleus gyromagnetic ratio

% Use Steck's table
gJ_s12=2.00233113;
gJ_p12=2/3;
gJ_p32=1.3362;

% Calculate from first principles
gJ_s12=gJ(0,1/2);
gJ_p12=gJ(1,1/2);
gJ_p32=gJ(1,3/2);

%% I MATRICES
% The nuclear spin I=3/2. Specifying raising and lowering operators in terms
% of elements to save space.

      
% Spin matrices for I=3/2 (COPY FOR J=3/2)
I_minus=  [0              0               0               0;
            low(3/2,3/2)  0               0               0;
            0               low(3/2,1/2)  0               0;
            0               0               low(3/2,-1/2) 0];
        
I_plus=   [0              up(3/2,1/2)  0               0;
             0              0               up(3/2,-1/2) 0;
             0              0               0               up(3/2,-3/2);
             0              0               0               0]; 
   
I_z=diag(3/2:-1:-3/2);
     
%% J MATRICES
  
% Spin matrices for J=1/2
J12_minus=[ 0               0;
            low(1/2,1/2)  0];
        
J12_plus=[  0               up(1/2,-1/2);
            0               0             ]; 
               
J12_z=diag(1/2:-1:-1/2);

      
% Spin matrices for J=3/2
J32_minus=  [0              0               0               0;
            low(3/2,3/2)  0               0               0;
            0               low(3/2,1/2)  0               0;
            0               0               low(3/2,-1/2) 0];
        
J32_plus=   [0              up(3/2,1/2)  0               0;
             0              0               up(3/2,-1/2) 0;
             0              0               0               up(3/2,-3/2);
             0              0               0               0]; 
   
J32_z=diag(3/2:-1:-3/2);


%% 4S_{1/2} Hamiltonian
JzIz=kron(J12_z,I_z);            % Jz*Iz
JpIm=kron(J12_plus,I_minus);     % J+*I- 
JmIp=kron(J12_minus,I_plus);     % J-*I+
JiIz=kron(eye(2*1/2+1),I_z);     % I_J*Iz
JzIi=kron(J12_z,eye(2*I+1));     % Jz*I_I   
JdotI=JzIz+0.5*(JpIm+JmIp);      % JdotI

Hs12=@(B) a_s12*JdotI+...
    muB*gI*JiIz*B+muB*gJ_s12*JzIi*B;


%% 4P_{1/2} Hamiltonian
JzIz=kron(J12_z,I_z);            % Jz*Iz
JpIm=kron(J12_plus,I_minus);     % J+*I- 
JmIp=kron(J12_minus,I_plus);     % J-*I+
JiIz=kron(eye(2*1/2+1),I_z);     % I_J*Iz
JzIi=kron(J12_z,eye(2*I+1));     % Jz*I_I   
JdotI=JzIz+0.5*(JpIm+JmIp);      % JdotI


Hp12=@(B) a_p12*JdotI+...
    muB*gI*JiIz*B+muB*gJ_p12*JzIi*B;

%% 4P_{3/2} Hamiltonian
JzIz=kron(J32_z,I_z);            % Jz*Iz
JpIm=kron(J32_plus,I_minus);     % J+*I- 
JmIp=kron(J32_minus,I_plus);     % J-*I+
JiIz=kron(eye(2*3/2+1),I_z);     % I_J*Iz
JzIi=kron(J32_z,eye(2*I+1));     % Jz*I_I   
JdotI=JzIz+0.5*(JpIm+JmIp);      % JdotI

% Jsquared
J2=J32_z^2+(0.5*(J32_plus+J32_minus))^2+(0.5/1i*(J32_plus-J32_minus))^2;
% Isquared
I2=I_z^2+(0.5*(I_plus+I_minus))^2+(0.5/1i*(I_plus-I_minus))^2;

% J^2I^2 
J2I2=kron(J2,I2); 

% Some vectors to test
vJ=[1; 0; 0; 0]; %mj=3/2;
vI=[1; 0; 0; 0; 0; 0; 0; 0; 0]; %mI=4;
vJI=kron(vJ,vI);% mJ=3/2;mI=4


J=3/2;
Hp32=@(B) a_p32*JdotI+...
    b_p32*(3*JdotI^2+1.5*JdotI-J2I2)/(2*I*(2*I-1)*J*(2*J-1))+...
    muB*gI*JiIz*B+muB*gJ_p32*JzIi*B;

%% Solve for eigenvalues and vectors

% Initialize energy matrix
Es12=zeros(size(Hs12(0),1),length(Bvec));
Ep12=zeros(size(Hp12(0),1),length(Bvec));
Ep32=zeros(size(Hp32(0),1),length(Bvec));

% Iterate over all magnetic fields
for kk=1:length(Bvec)
    % Find the hamiltonians
    Hs12_mat=Hs12(Bvec(kk));
    Hp12_mat=Hp12(Bvec(kk));
    Hp32_mat=Hp32(Bvec(kk));
    
    % Solve the hamiltonians
    Es12(:,kk)=sort(eig(Hs12_mat));
    Ep12(:,kk)=sort(eig(Hp12_mat));
    Ep32(:,kk)=sort(eig(Hp32_mat));
end

% Assign output data
out.Es12=Es12;
out.Ep12=Ep12;
out.Ep32=Ep32;

%% Plot Results
% Plot the results.  Limits are chosen to match Tiecke. Also adding labels
% for F, mJ and the zero field splittings.
% 5S_1/2
hF1=figure(2000);
clf
set(hF1,'color','w','Name','Rb 5 S 1/2');
hold on
for kk=1:size(Es12,1)
   plot(Bvec,Es12(kk,:)*1e-3,'k-','linewidth',1) 
end
xlim([0 15000]);
ylim([-25 25]);

hF1.Position(3:4)=[600 300];
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (GHz)');

text(0.02,.98,'$^{87}\mathrm{Rb}~5\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','top','fontsize',18);

text(1.005,.9,'$m_\mathrm{J}=+\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(1.005,.05,'$m_\mathrm{J}=-\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(.15,.7,'$F=2$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(.15,.3,'$F=1$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(0,1E-3*Es12(1,1),[num2str(Es12(1,1)*1E-3,'%.3f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(0,1E-3*Es12(end,1),[num2str(Es12(end,1)*1E-3,'%.3f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');

% 5P_1/2
hF2=figure(2001);
clf
set(hF2,'color','w','name','Rb 5 P 1/2');
hold on
for kk=1:size(Ep12,1)
   plot(Bvec,Ep12(kk,:),'k-','linewidth',1) 
end
xlim([0 5000]);
ylim([-2500 2500]);
hF2.Position(3:4)=[600 300];
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
text(0.02,.98,'$^{87}\mathrm{Rb}~5\mathrm{P}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','top','fontsize',18);
% 
text(1.005,.95,'$m_\mathrm{J}=+\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(1.005,.05,'$m_\mathrm{J}=-\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(.15,.7,'$F=2$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(.15,.25,'$F=1$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(0,Ep12(1,1),[num2str(Ep12(1,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(0,Ep12(end,1)-15,[num2str(Ep12(end,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');

% 5P_3/2
hF3=figure(2002);
clf
set(hF3,'color','w','name','87Rb P 3/2');
hold on
for kk=1:size(Ep32,1)
   plot(Bvec,Ep32(kk,:),'k-','linewidth',1) 
end
xlim([0 500]);
ylim([-1500 1500]);
hF3.Position(3:4)=[600 300];
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
text(0.02,.98,'$^{87}\mathrm{Rb}~5\mathrm{P}_{3/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','top','fontsize',18);

text(.85,.75,'$m_\mathrm{J}=+\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(.85,.25,'$m_\mathrm{J}=-\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(.4,.85,'$m_\mathrm{J}=+\frac{3}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(.4,.15,'$m_\mathrm{J}=-\frac{3}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(.13,.7,'$F=3$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(.13,.3,'$F=0$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(-6*0,Ep32(1,1)-30,[num2str(Ep32(1,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(-6*0,Ep32(2,1)+20,[num2str(Ep32(2,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(-6*0,Ep32(5,1)-20,[num2str(Ep32(5,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(0,Ep32(end,1)+20,[num2str(Ep32(end,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
%% Helper functions

    % Coeffcients for raising operator     
    function coeff=up(jj,m)
        coeff=sqrt(jj*(jj+1)-m*(m+1));
    end

    % Coefficients for lowering operator
    function coeff=low(jj,m)
        coeff=sqrt(jj*(jj+1)-m*(m-1));
    end

    function out=gJ(L,J)
    %    out=3/2+(S*(S+1)-L*(L+1)) /(2*J*(J+1));
       %https://en.wikipedia.org/wiki/Land%C3%A9_g-factor
        out=gL*(J*(J+1)-S*(S+1)+L*(L+1))/(2*J*(J+1))+...
        gS*(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1));

    end
end


