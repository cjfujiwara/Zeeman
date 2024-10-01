function [out,hF1,hF2,hF3]=K_zeeman
% Author : C Fujiwara
%
% Calculate Zeeman shifts for 40K
%
% Output the eigen energies for a variety of magnetic fields Bvec. The
% output structure "out" is assigned the eigen values
out=struct;
Bvec=linspace(0,1000,1E4);
% Bvec = [Bvec 1e5];
Bvec(1)=[];
out.B=Bvec;


%% Constants
% Bohr magneton in MHz/Gauss
muB=1.39962449 ;    

% 40K has S=1/2 and I=4;
S=0.5;
I=4;

%% Coupling Constants
a_s12=-285.7308;    % 4S_1/2 magnetic dipole hyperfine
a_p12=-34.523;      % 4P_1/2 magnetic dipole hyperfine
a_p32=-7.585;       % 4P_3/2 magnetic dipole hyperfine
b_p32=-3.445;       % 4P_3/2 electric quadrupole hyperfine

gS=2.0023193043622; % Electron gyromagnetic ratio
gL=1;               % Electron orbital gyromagnetic ratio
gI=0.000176490;     % 40K nucleus gyromagnetic ratio

% Use Tiecke's table
gJ_s12=2.00229421;
gJ_p12=2/3;
gJ_p32=4/3;

% Calculate from first principles (gets numbers consistent with level
% diagrams)
gJ_s12=gJ(0,1/2);
gJ_p12=gJ(1,1/2);
gJ_p32=gJ(1,3/2);


%% I MATRICES
% The nuclear spin I=4. Specifying raising and lowering operators in terms
% of elements to save space.

% Lowering operator
I_minus=zeros(2*I+1,2*I+1);
I_minus(2,1)=low(I,4);
I_minus(3,2)=low(I,3);
I_minus(4,3)=low(I,2);
I_minus(5,4)=low(I,1);
I_minus(6,5)=low(I,0);
I_minus(7,6)=low(I,-1);
I_minus(8,7)=low(I,-2);
I_minus(9,8)=low(I,-3);

% Raising Operater
I_plus=zeros(2*I+1,2*I+1);
I_plus(1,2)=up(I,3);
I_plus(2,3)=up(I,2);
I_plus(3,4)=up(I,1);
I_plus(4,5)=up(I,0);
I_plus(5,6)=up(I,-1);
I_plus(6,7)=up(I,-2);
I_plus(7,8)=up(I,-3);
I_plus(8,9)=up(I,-4);

% Z 
I_z=diag(4:-1:-4);


     
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
JzIz_S12=kron(J12_z,I_z);                       % Jz*Iz
JpIm_S12=kron(J12_plus,I_minus);                % J+*I- 
JmIp_S12=kron(J12_minus,I_plus);                % J-*I+
JiIz_S12=kron(eye(2*1/2+1),I_z);                % I_J*Iz
JzIi_S12=kron(J12_z,eye(2*I+1));                % Jz*I_I   
JdotI_S12=JzIz_S12+0.5*(JpIm_S12+JmIp_S12);      % JdotI

Hs12=@(B) a_s12*JdotI_S12+...
    muB*gI*JiIz_S12*B+muB*gJ_s12*JzIi_S12*B;

%% 4P_{1/2} Hamiltonian
JzIz_P12=kron(J12_z,I_z);            % Jz*Iz
JpIm_P12=kron(J12_plus,I_minus);     % J+*I- 
JmIp_P12=kron(J12_minus,I_plus);     % J-*I+
JiIz_P12=kron(eye(2*1/2+1),I_z);     % I_J*Iz
JzIi_P12=kron(J12_z,eye(2*I+1));     % Jz*I_I   
JdotI_P12=JzIz_P12+0.5*(JpIm_P12+JmIp_P12);      % JdotI


Hp12=@(B) a_p12*JdotI_P12+...
    muB*gI*JiIz_P12*B+muB*gJ_p12*JzIi_P12*B;

% Some vectors to test
vJ=[1; 0]; %mj=1/2;
vI=[1; 0; 0; 0; 0; 0; 0; 0; 0]; %mI=4;
vJI_P12=kron(vJ,vI);% mJ=1/2;mI=4

%% 4P_{3/2} Hamiltonian
JzIz_P32=kron(J32_z,I_z);            % Jz*Iz
JpIm_P32=kron(J32_plus,I_minus);     % J+*I- 
JmIp_P32=kron(J32_minus,I_plus);     % J-*I+
JiIz_P32=kron(eye(2*3/2+1),I_z);     % I_J*Iz
JzIi_P32=kron(J32_z,eye(2*I+1));     % Jz*I_I   
JdotI_P32=JzIz_P32+0.5*(JpIm_P32+JmIp_P32);      % JdotI

% Jsquared
J2=J32_z^2+(0.5*(J32_plus+J32_minus))^2+(0.5/1i*(J32_plus-J32_minus))^2;
% Isquared
I2=I_z^2+(0.5*(I_plus+I_minus))^2+(0.5/1i*(I_plus-I_minus))^2;

% J^2I^2 
J2I2=kron(J2,I2); 

% Some vectors to test
vJ=[1; 0; 0; 0]; %mj=3/2;
vI=[1; 0; 0; 0; 0; 0; 0; 0; 0]; %mI=4;
vJI_P32=kron(vJ,vI);% mJ=3/2;mI=4


J=3/2;
Hp32=@(B) a_p32*JdotI_P32+...
    b_p32*(3*JdotI_P32^2+1.5*JdotI_P32-J2I2)/(2*I*(2*I-1)*J*(2*J-1))+...
    muB*gI*JiIz_P32*B+muB*gJ_p32*JzIi_P32*B;

%% Solve for eigenvalues and vectors

% Initialize energy matrix
Es12=zeros(size(Hs12(0),1),length(Bvec));
Ep12=zeros(size(Hp12(0),1),length(Bvec));
Ep32=zeros(size(Hp32(0),1),length(Bvec));

vs12=zeros(size(Hs12(0),1),size(Hs12(0),1),length(Bvec));
vp12=zeros(size(Hp12(0),1),size(Hp12(0),1),length(Bvec));
vp32=zeros(size(Hp32(0),1),size(Hp32(0),1),length(Bvec));

% Iterate over all magnetic fields
for kk=1:length(Bvec)
    % Find the hamiltonians

    Hs12_mat=Hs12(Bvec(kk));
    Hp12_mat=Hp12(Bvec(kk));
    Hp32_mat=Hp32(Bvec(kk));
    
    % 4S 1/2
    [v,d]=eig(Hs12_mat);
     [d,inds]=sort(diag(d));v=v(:,inds);        
    Es12(:,kk)=d;   

    vs12(:,:,kk)=v;      
    
    
    
    % 4P 1/2
    [v,d]=eig(Hp12_mat);
    [d,inds]=sort(diag(d));v=v(:,inds);    
    Ep12(:,kk)=d;
    vp12(:,:,kk)=v;    
    
    % 4P 3/2
    [v,d]=eig(Hp32_mat);
    [d,inds]=sort(diag(d));v=v(:,inds);    
    Ep32(:,kk)=d;
    vp32(:,:,kk)=v;
    
  
end

% Calculate mI and mJ
mImJ_S12 = zeros(size(vs12,1),2);
for kk = 1:size(vs12,2)
    mI = vs12(:,kk,end)'*JiIz_S12*vs12(:,kk,end);
    mJ = vs12(:,kk,end)'*JzIi_S12*vs12(:,kk,end); 
    mImJ_S12(kk,:) = [mI mJ];
end

mImJ_P12 = zeros(size(vp12,1),2);
for kk = 1:size(vp12,2)
    mI = vp12(:,kk,end)'*JiIz_P12*vp12(:,kk,end);
    mJ = vp12(:,kk,end)'*JzIi_P12*vp12(:,kk,end); 
    mImJ_P12(kk,:) = [mI mJ];
end

mImJ_P32 = zeros(size(vp32,1),2);
for kk = 1:size(vp32,2)
    mI = vp32(:,kk,end)'*JiIz_P32*vp32(:,kk,end);
    mJ = vp32(:,kk,end)'*JzIi_P32*vp32(:,kk,end); 
    mImJ_P32(kk,:) = [mI mJ];
end

% Calculate Basis
basis_mImJ_S12 = zeros(size(vs12,1),2);
for kk = 1:size(vs12,2)
    v=zeros(size(vs12,1),1);
    v(kk)=1;
    mI = v'*JiIz_S12*v;
    mJ = v'*JzIi_S12*v; 
    basis_mImJ_S12(kk,:) = [mI mJ];
end

% Calculate Basis
basis_mImJ_P12 = zeros(size(vp12,1),2);
for kk = 1:size(vp12,2)
    v=zeros(size(vp12,1),1);
    v(kk)=1;
    mI = v'*JiIz_P12*v;
    mJ = v'*JzIi_P12*v; 
    basis_mImJ_P12(kk,:) = [mI mJ];
end


% Calculate Basis
basis_mImJ_P32 = zeros(size(vp32,1),2);
for kk = 1:size(vp32,2)
    v=zeros(size(vp32,1),1);
    v(kk)=1;
    mI = v'*JiIz_P32*v;
    mJ = v'*JzIi_P32*v; 
    basis_mImJ_P32(kk,:) = [mI mJ];
end

% Assign output data
out.Es12=Es12;
out.vs12=vs12;
out.mImJ_S12 = mImJ_S12;
out.basis_mImJ_S12 = basis_mImJ_S12;

% out.JiIz_S12 = JiIz_S12; % Iz operator for 4S_1/2
% out.JzIi_S12= JzIi_S12;  % Jz operator for 4S_1/2

out.Ep12=Ep12;
out.vp12=vp12;
out.mImJ_P12 = mImJ_P12;
out.basis_mImJ_P12 = basis_mImJ_P12;

% out.JiIz_P12 = JiIz_P12; % Iz operator for 4P_1/2
% out.JzIi_P12 = JzIi_P12;  % Jz operator for 4P_1/2

out.Ep32=Ep32;
out.vp32=vp32;
out.mImJ_P32 = mImJ_P32;

out.basis_mImJ_P32 = basis_mImJ_P32;

% 
% out.JiIz_P32 = JiIz_P32; % Iz operator for 4P_3/2
% out.JzIi_P32 = JzIi_P32;  % Jz operator for 4P_3/2

%% Plot Results
% Plot the results.  Limits are chosen to match Tiecke. Also adding labels
% for F, mJ and the zero field splittings.
% 4S_1/2
hF1=figure(1000);
clf
set(hF1,'color','w');
hold on
for kk=1:size(Es12,1)
   plot(Bvec,Es12(kk,:),'k-','linewidth',1) 
end
xlim([0 1000]);
ylim([-2000 2000]);
hF1.Position(3:4)=[600 400];
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
text(0.02,.98,'$^{40}\mathrm{K}~4\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','top','fontsize',18);

text(1.005,.85,'$m_\mathrm{J}=+\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(1.005,.15,'$m_\mathrm{J}=-\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(.15,.8,'$F=7/2$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(.15,.2,'$F=9/2$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(0,Es12(1,1),[num2str(Es12(1,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(0,Es12(end,1),[num2str(Es12(end,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');


text(.78,.04,['$m_I = ' num2str(mImJ_S12(1,1)) '$'],'fontsize',8,'units','normalized',...
    'horizontalalignment','right','interpreter','latex');
% text(.78,.4,['$m_I = ' num2str(mImJ_S12(9,1)) '$'],'fontsize',8,'units','normalized',...
%     'horizontalalignment','right','interpreter','latex');


text(.78,.58,['$m_I = ' num2str(mImJ_S12(10,1)) '$'],'fontsize',8,'units','normalized',...
    'horizontalalignment','right','interpreter','latex');
% text(.78,.97,['$m_I = ' num2str(mImJ_S12(end,1)) '$'],'fontsize',8,'units','normalized',...
%     'horizontalalignment','right','interpreter','latex');
% 4P_1/2
hF2=figure(1001);
clf
set(hF2,'color','w');
hold on
for kk=1:size(Ep12,1)
   plot(Bvec,Ep12(kk,:),'k-','linewidth',1) 
end
xlim([0 400]);
ylim([-300 300]);
hF2.Position(3:4)=[600 400];
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
text(0.02,.98,'$^{40}\mathrm{K}~4\mathrm{P}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','top','fontsize',18);

text(1.005,.85,'$m_\mathrm{J}=+\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(1.005,.15,'$m_\mathrm{J}=-\frac{1}{2}$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(.15,.8,'$F=7/2$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(.15,.2,'$F=9/2$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(0,Ep12(1,1),[num2str(Ep12(1,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(0,Ep12(end,1)-15,[num2str(Ep12(end,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');


text(.78,.08,['$m_I = ' num2str(mImJ_P12(1,1)) '$'],'fontsize',8,'units','normalized',...
    'horizontalalignment','right','interpreter','latex');
% text(.78,.4,['$m_I = ' num2str(mImJ_P12(9,1)) '$'],'fontsize',8,'units','normalized',...
%     'horizontalalignment','right','interpreter','latex');


text(.78,.58,['$m_I = ' num2str(mImJ_P12(10,1)) '$'],'fontsize',8,'units','normalized',...
    'horizontalalignment','right','interpreter','latex');

% text(.78,.9,['$m_I = ' num2str(mImJ_P12(end,1)) '$'],'fontsize',8,'units','normalized',...
%     'horizontalalignment','right','interpreter','latex');


% 4P_3/2
hF3=figure(1002);
clf
set(hF3,'color','w');
hold on
for kk=1:size(Ep32,1)
   plot(Bvec,Ep32(kk,:),'k-','linewidth',1) 
end
xlim([0 200]);
ylim([-400 400]);
hF3.Position(3:4)=[600 400];
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
text(0.02,.98,'$^{40}\mathrm{K}~4\mathrm{P}_{3/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','top','fontsize',18);

text(.85,.65,'$m_\mathrm{J}=+\frac{1}{2}$','interpreter','latex','fontsize',12,'units','normalized',...
    'horizontalalignment','left');
text(.85,.35,'$m_\mathrm{J}=-\frac{1}{2}$','interpreter','latex','fontsize',12,'units','normalized',...
    'horizontalalignment','left');

text(.55,.8,'$m_\mathrm{J}=+\frac{3}{2}$','interpreter','latex','fontsize',12,'units','normalized',...
    'horizontalalignment','left');
text(.55,.2,'$m_\mathrm{J}=-\frac{3}{2}$','interpreter','latex','fontsize',12,'units','normalized',...
    'horizontalalignment','left');

text(.15,.8,'$F=5/2$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');
text(.15,.2,'$F=11/2$','interpreter','latex','fontsize',10,'units','normalized',...
    'horizontalalignment','left');

text(-6,Ep32(1,1),[num2str(Ep32(1,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(-6,Ep32(13,1),[num2str(Ep32(13,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(-6,Ep32(end-7,1),[num2str(Ep32(end-7,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');
text(-6,Ep32(end,1)+10,[num2str(Ep32(end,1),'%.1f') ' '],'fontsize',7,'units','data',...
    'horizontalalignment','right');



text(.35,.15,['$' num2str(mImJ_P32(1,1)) '$'],'fontsize',8,'units','normalized',...
    'horizontalalignment','right','interpreter','latex');
text(.35,.37,['$' num2str(mImJ_P32(10,1)) '$'],'fontsize',8,'units','normalized',...
    'horizontalalignment','right','interpreter','latex');
text(.35,.52,['$' num2str(mImJ_P32(19,1)) '$'],'fontsize',8,'units','normalized',...
    'horizontalalignment','right','interpreter','latex');
text(.35,.65,['$' num2str(mImJ_P32(19,1)) '$'],'fontsize',8,'units','normalized',...
    'horizontalalignment','right','interpreter','latex');
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


