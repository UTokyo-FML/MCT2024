function [Lneg,Fb,gelF,Qk,udot1] = variableHyperGel2_variableDt_passive(TCatilde,TCastar,Tstar,X_w,X_p, dt,Qk,L,CL)
% gelF corresponds to the force/area exerted by gel on the cell; equals to
% <gelF> = <"cross-bridge force" - "passive force">
% units: s, um, nN, kPa ( nN / um^2 = kPa )
Ap=1008e+04; %[kPa/um/mM]
Aw=Ap/5;  %[kPa/um/mM]
% Maxwell Model
lgel = 100; % [um] length of the gel
Sgel = 600;% [um^2] Myocyte cut area is ~18*18. I assume ~2 times of that is affected.
lcell= 65; % [um] Cell length. Symmetricity reduces the analysis field to 1/2.
Scell= 324;% [um^2] Cell size in the experiment was 131.2 * 18.4.
gammaNB = 0.6; % ratio
mxiter  = 50;
L0=0.97;
Nmaxwell = 13; % Same as N_tau. Number of maxwell elements
[Egel,eta] = setparameter(Nmaxwell,CL); % E [kPa], eta [s]
[Bgel,Bk,betak]= calcMaxwellStiffMat_1dim(Egel,eta,dt,lgel,Sgel,Nmaxwell);

udot0 = 0;
u0=(L/L0-1.0)*lcell;
%L=(u0/lcell+1)*L0
% iteration 0 update
udot1=udot0;
u1 = u0 + dt*udot0;
Lneg = (u1/lcell+1)*L0;
Fb =Aw*TCatilde*(Lneg-X_w) + Ap*(TCastar+Tstar)*(Lneg-X_p);

alpha = u1/lcell+1.0;
% Rin-Yin Model. The function is nonlinear thus called at each iteration
[S,Dcell]= HumphreyPassive(alpha);
[Kcell,equivf] = element1D(S,Dcell,Scell,lcell,alpha);
Pmxwl = 2*Bgel*udot0+betak*Qk;
for iter = 1:mxiter % Newton Raphson iteration
    % Merge, extracting node 2 equation only (node 1 and 3 are bounded)
    resid = -Fb*Scell-(-equivf)-Pmxwl;
    if(resid^2 < 1.E-16)
        break;
    end
    astiff = gammaNB*dt*(Kcell+(Aw*TCatilde+ Ap*(TCastar+Tstar))*Scell)+ Bgel;
    dudot = resid/astiff;
    udot1 = udot1 + dudot;
    u1    = u1 + gammaNB*dt*dudot;

    Pmxwl = Bgel*(udot1+udot0)+betak*Qk;
    alpha = u1/lcell+1.0;
    Lneg = alpha*L0;
    Fb =Aw*TCatilde*(Lneg-X_w) + Ap*(TCastar+Tstar)*(Lneg-X_p);
    [S,Dcell]= HumphreyPassive(alpha);
    [Kcell,equivf] = element1D(S,Dcell,Scell,lcell,alpha);
end
Qk = Bk*(udot1+udot0)+ betak.' .*Qk; % for each k. don't take sum
gelF=-Pmxwl/Scell;

function initGel(dt)
% Those parameters never changes if dt is constant. The function is called once.
Nmaxwell = 13; % Same as N_tau. Number of maxwell elements
[Egel,eta] = setparameter(Nmaxwell); % E [kPa], eta [s]
[Bgel,Bk,betak]= calcMaxwellStiffMat_1dim(Egel,eta,dt,lgel,Sgel,Nmaxwell);

function [E,eta] = setparameter(Nmaxwell,CL)
% Mohammad paper viscoelastic parameters (Eshelby Inclusion paper)
zeta0 = 3.36224; 
g0   = 14.7; %[kPa]
g1   = 34.4;
g2   = 0.0866;
G0   = g0/2*(1+tanh(g1*(CL-g2))); %eq 9a
h0   = 0.00486; %[s]
h1   = 1.17; %[s]
h2   = -2.62; %[s]
tau0 = h0 + h1*CL+h2*CL^2; %eq 9b
s0   = 0.58;
s1   = 55;
sigma0 = s0*(1-exp(-s1*CL)); %eq 9c
% G0
% tau0
% sigma0
% pause(10)
Ginfty= 0;     % [kPa] G0*0.1; %
zeta = zeros(1,Nmaxwell);
tau = zeros(1,Nmaxwell);
E   = zeros(1,Nmaxwell);
eta = zeros(1,Nmaxwell);
for i = 1:Nmaxwell
    tau(i) = 2^(i+2-Nmaxwell);
    zeta(i)= zeta0*exp(-(log(tau(i)/tau0)/sigma0)^2);
    E(i)   = 3*(G0-Ginfty)*zeta(i);
    eta(i) = E(i)*tau(i);
end
% for debugging, check G(t): relaxation function
% Mohammad et al. Int. J. Eng. Sci 165 (2021) figure 2, (b) and (c)
% Amplitude is slightly different though....
tsteps = -2: 0.05 :2;
omega  = zeros(size(tsteps));
freq   = zeros(size(tsteps));
Gstorage = zeros(size(tsteps));
Gloss    = zeros(size(tsteps));
for i= 1:size(tsteps,2)
    freq(i)  = 10^tsteps(i);
    omega(i) = 2*pi*freq(i);
    Gstorage(i) = Ginfty;
    Gloss(i)    = 0;
    for j=1:Nmaxwell
        omegatau  = omega(i)*tau(j);
        omegatau2 = omegatau^2;
        denom     = 1+omegatau2;
        Gstorage(i) = Gstorage(i)+ (G0-Ginfty)*omegatau2*zeta(j)/denom;
        Gloss(i)    = Gloss(i)   + (G0-Ginfty)*omegatau *zeta(j)/denom;
    end
end
% figure(1);
% subplot(2,1,1); semilogx(freq,Gstorage); title('G storage');
% subplot(2,1,2); semilogx(freq,Gloss); title('G loss');


function [B,Bk,betak] = calcMaxwellStiffMat_1dim(E,eta,dt,d,S,Nmaxwell)
% Model size affects the stiffness matrix. 
% d: length of the FE element
% S: cut area of the FE element
zeta = S/d*eta; %Corresponding viscosity coefficient
gamma= S/d*E;   %Corresponding spring constant
B=0;
Bk = zeros(Nmaxwell,1);
betak = zeros(1,Nmaxwell);
for i = 1:Nmaxwell
    denom = 2/(gamma(i)*dt)+1/zeta(i);
    Bk(i) = 1/denom;
    B = B + Bk(i);
    betak(i) = (2/(gamma(i)*dt)-1/zeta(i)) /denom;
end
if(eta == 0 & E == 0)
    betak =zeros(1,Nmaxwell);
    Bk = zeros(Nmaxwell,1);
end

function [Kcell,equivf] = element1D(S,D,Area,lcell,alpha)
% S: 2nd Piola Kirchhoff Stress (component 1,1)
% D: = 4 ddW/dCdC
% Area: Element cut area
% alpha: stretch (1+epsilon) to the fiber direction (1D length direction)
%
% The function culculates the following K matrix (1,1) component value
%  {del u} (D*alapha^2/L^2 + S/L^2) |  1  -1 | Vol {delta u}
%                                   | -1   1 |
% Vol = Scell * L
%  {del u} (D*alapha^2 + S)/L |  1  -1 | Area {delta u}  = | K -K |
%                             | -1   1 |                   |-K  K |
% equivalent force
% [B]^T [S]*Vol =  alpha/L * [-1, 1]*S*Vol = alpha*Area*S [-1, 1]

Kcell = (alpha^2*D+S)*Area/lcell;
equivf = -alpha*Area*S;

function [Kcell,equivf] = NegroniPassive(alpha,Scell,Lcell)
Ke=105000; %[kPa/um^5]
Le=10;
L0=0.97;
x = L0*alpha;
Kcell = 5*Ke*(x-L0)^4+Le;
Kcell = Kcell*Scell/Lcell;
equivf = Ke*(x-L0)^5+Le*(x-L0);
equivf = -equivf*Scell/Lcell;

function [Kcell,equivf] = GelPassiveSpring(u,S,L,CL)
Ke=1050000*CL; %[kPa/um^5]
Le=10*CL;
dudx = u/L;
Kcell = 5*Ke*dudx^4+Le;
Kcell = Kcell*S;
equivf = Ke*dudx^5+Le*dudx;
equivf = -equivf*S;

function [S,Dmat] = HumphreyPassive(alpha)
%unit kPa, um, nN
cp1 =  3080.D-3;
cp2 =  3240.D-3;
cp3 =   359.D-3;
cp4 = -1940.D-3;
cp5 =  1660.D-3;

% Assume 1D stretch alpha, incompressible. Then C mat becomes
%     | alpha^2                   |
% C = |          1/alpha          |
%     |                   1/alpha |
pi1 = alpha^2+2/alpha;
pi4 = alpha^2;
ri4 = pi4;
ri1 = alpha^2+2/alpha;
ri1m = ri1 -3.D0;
ri4m = alpha - 1.D0;
ri4ma = ri4 - 1.D0;

% passive
dwdr1 = cp3 + cp4*ri4m + 2.D0*cp5*ri1m;
dwdr4 = (1.D0/2.D0)/alpha *(2.D0*cp1*ri4m +3.D0*cp2*(ri4m)^ ...
    2+cp4*ri1m);
dwd11 = 2.D0*cp5;
dwd14 = (1.D0/2.D0)/alpha*cp4;
dwd41 = (1.D0/2.D0)/alpha*cp4;
dwd44 = (-1.D0/4.D0)*((ri4)^(-3.D0/2.D0)) *(2.D0*cp1*ri4m+3.D0*cp2*(ri4m)^2+cp4*ri1m) ...
       + (1.D0/2.D0)/alpha*(cp1/alpha+3.D0*cp2/alpha*ri4m);

% i,j,k,l = 1
dr1dcij = 2/3*(1-1/alpha^2); % = 1.D0/rdetc1 * (del(i,j)-pi1*rci(i,j)/3.D0);
dr4dcij = 1; % = dn(i)*dn(j)
S       = 2.D0 * ( dwdr1 * dr1dcij + dwdr4 * dr4dcij );
dr1dckl = 2/3*(1-1/alpha^2); % = 1.D0/rdetc1 * (del(k,l)-pi1*rci(k,l)/3.D0)
dr4dckl = 1; % = dn(k)*dn(l)
ddr1dcdc = -1/3*(1-pi1/alpha^2/3)/alpha^2+ (-1/alpha^2+pi1/alpha^2/alpha^2)/3.0; % = -1.D0/rdetc1/3.D0 * (del(i,j)-pi1*rci(i,j)/3.D0)*rci(k,l) + 1.D0/rdetc1 * (-rci(i,j)*del(k,l) + pi1*rci(i,k)*rci(l,j)) /3.D0
ddr4dcdc = 0.D0;
Dmat = 4.0 *((dwd11*dr1dckl + dwd14*dr4dckl) * dr1dcij  + dwdr1 * ddr1dcdc ...
            +(dwd41*dr1dckl + dwd44*dr4dckl) * dr4dcij  + dwdr4 * ddr4dcdc );

function output = DebugHumphreyCheck()
% Debug check the stress-stretch relationship
alpha = 0.6 : 0.01: 2;
fr = 0: 0.2: 1;

Sarray = zeros(size(alpha,2));
for i=1:size(alpha,2)
   [Sarray(i),Dmat] = HumphreyPassive(alpha(i));
end
figure(2);
%plot(alpha,Sarray(:,1),"-",alpha,Sarray(:,3),"--",alpha,Sarray(:,5),"-.",alpha,Sarray(:,6),":");
plot(alpha,Sarray);
xlabel('strain')
ylabel('stress')
output = 1;

function [S,Dmat] = RinYin(fr, alpha)
%unit kPa, um, nN
cp1 =  3080.D-3;
cp2 =  3240.D-3;
cp3 =   359.D-3;
cp4 = -1940.D-3;
cp5 =  1660.D-3;
ca0 =  0.D-3 * fr;
ca1 = -0.638903D1 * fr;
ca2 =  0.179707D2 * fr;
ca3 =  0.173676D2 * fr;
ca4 =  0.760996D1 * fr^2;
%    ca5 =  0.194223D2 * fr**2;
ca5 =  0.D0;

% Assume 1D stretch alpha, incompressible. Then C mat becomes
%     | alpha^2                   |
% C = |          1/alpha          |
%     |                   1/alpha |
pi1 = alpha^2+2/alpha;
pi4 = alpha^2;
ri4 = pi4;
ri1 = alpha^2+2/alpha;
ri1m = ri1 -3.D0;
ri4m = alpha - 1.D0;
ri4ma = ri4 - 1.D0;

% passive
dwdr1p = cp3 + cp4*ri4m + 2.D0*cp5*ri1m;
dwdr4p = (1.D0/2.D0)/alpha *(2.D0*cp1*ri4m +3.D0*cp2*(ri4m)^ ...
    2+cp4*ri1m);
dwd11p = 2.D0*cp5;
dwd14p = (1.D0/2.D0)/alpha*cp4;
dwd41p = (1.D0/2.D0)/alpha*cp4;
dwd44p = (-1.D0/4.D0)*((ri4)^(-3.D0/2.D0)) *(2.D0*cp1*ri4m+3.D0*cp2*(ri4m)^2+cp4*ri1m) ...
       + (1.D0/2.D0)/alpha*(cp1/alpha+3.D0*cp2/alpha*ri4m);

% active
dwdr1a = ca1 * ri4ma + 2.D0 * ca2 * ri1m + ca4;
dwdr4a = ca1 * ri1m + 2.D0 * ca3 * ri4ma + ca5;
dwd11a = 2.D0 * ca2;
dwd14a = 1.D0 * ca1;
dwd41a = 1.D0 * ca1;
dwd44a = 2.D0 * ca3;

dwdr1 = dwdr1p + dwdr1a;
dwdr4 = dwdr4p + dwdr4a;
dwd11 = dwd11p + dwd11a;
dwd14 = dwd14p + dwd14a;
dwd41 = dwd41p + dwd41a;
dwd44 = dwd44p + dwd44a;

% i,j,k,l = 1
dr1dcij = 2/3*(1-1/alpha^2); % = 1.D0/rdetc1 * (del(i,j)-pi1*rci(i,j)/3.D0);
dr4dcij = 1; % = dn(i)*dn(j)
S       = 2.D0 * ( dwdr1 * dr1dcij + dwdr4 * dr4dcij );
dr1dckl = 2/3*(1-1/alpha^2); % = 1.D0/rdetc1 * (del(k,l)-pi1*rci(k,l)/3.D0)
dr4dckl = 1; % = dn(k)*dn(l)
ddr1dcdc = -1/3*(1-pi1/alpha^2/3)/alpha^2+ (-1/alpha^2+pi1/alpha^2/alpha^2)/3.0; % = -1.D0/rdetc1/3.D0 * (del(i,j)-pi1*rci(i,j)/3.D0)*rci(k,l) + 1.D0/rdetc1 * (-rci(i,j)*del(k,l) + pi1*rci(i,k)*rci(l,j)) /3.D0
ddr4dcdc = 0.D0;
Dmat = 4.0 *((dwd11*dr1dckl + dwd14*dr4dckl) * dr1dcij  + dwdr1 * ddr1dcdc ...
            +(dwd41*dr1dckl + dwd44*dr4dckl) * dr4dcij  + dwdr4 * ddr4dcdc );

function output = DebugRinYinCheck()
% Debug check the stress-stretch relationship
alpha = 0.6 : 0.01: 2;
fr = 0: 0.2: 1;

Sarray = zeros(size(alpha,2),size(fr,2));
for i=1:size(alpha,2)
    for j=1:size(fr,2)
        [Sarray(i,j),Dmat] = RinYin(fr(j), alpha(i));
    end
end
figure(2);
%plot(alpha,Sarray(:,1),"-",alpha,Sarray(:,3),"--",alpha,Sarray(:,5),"-.",alpha,Sarray(:,6),":");
plot(alpha,Sarray(:,1),alpha,Sarray(:,2),alpha,Sarray(:,3),alpha,Sarray(:,4),alpha,Sarray(:,5),alpha,Sarray(:,6));
legend("passive","fr = 0.2","fr = 0.4","fr = 0.6","fr = 0.8","fr = 1.0");
output = 1;