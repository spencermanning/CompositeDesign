%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spencer Manning
% ME 456 CLT CODE
% 4/12/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code does all of the following actions: 

% Accepts user inputs for:
%   No. plies, thicknesses, orientations
% Calculates:
%   Matrix and resin and volume fraction
%   Forces/moments

% Calculates Engineering parameters:
%   E, G, v
%   SL+,SL-, etc. 
%   Q
%   Qbar
%   ABBD Matrix

% Calculates responses:
%   Strain and curvature on middle plane
%   Strain on top and bottom of each ply in x-y frame
%   Stress on top and bottom of each ply in x-y frame
%   Stresses and strains on top and bottom of each ply in 1-2 frame

% Determines the failure mode:
%   Apply relevant failure model at top and bottom of each ply, and test for failure

clear
clc
close all

% Define the names of all the materials
NameF = {'Boron', 'HMS', 'AS', 't300', 'Kevlar', 'SG', 'EG'};
NameM = {'LM','IMLS','IMHS','HM','Polyimide','PMR'};

%%
% -----Fibers-----
%             fd(m)     pf(g/m^3)   E1f(GPa)  E2f(GPa)  G1f(GPa)  G2f(GPa)  v1f  v2f     SLpf(MPa)  SLmf(MPa) ef1
Boron =     [0.00014224   2629591    400       400       167       167       .2   .2      4140       4830      .008];  %Boron
HMS =       [7.62e-6      1937593    379       6.21      7.58      4.83      .2   .00327  1720       1380      .007];  %HMS
AS =        [7.62e-6      1743834    213.7     13.8      13.8      6.89      .2   .0129   2410       1790      .018];  %AS
t300 =      [7.62e-6      1771514    221       13.8      8.96      4.83      .2   .0125   2410       2070      .014];  %t300
Kevlar =    [1.1684e-5    1467035    152       4.14      2.90      1.52      .35  .00954  2760        517      .024];  %Kevlar
SG =        [9.144e-6     2491191    85.5       855       35.6      35.6      .2   .2      4140       0         .057];  %S-G
EG =        [9.144e-6     2491191    731       731       30.1      35.6      .22  .22     2760       0         .048];  %E-G

% Organize the material properties
fd =    [Boron(1) HMS(1) AS(1) t300(1) Kevlar(1) SG(1) EG(1)];
pf =    [Boron(2) HMS(2) AS(2) t300(2) Kevlar(2) SG(2) EG(2)];
E1f =   [Boron(3) HMS(3) AS(3) t300(3) Kevlar(3) SG(3) EG(3)];
E2f =   [Boron(4) HMS(4) AS(4) t300(4) Kevlar(4) SG(4) EG(4)];
G1f =   [Boron(5) HMS(5) AS(5) t300(5) Kevlar(5) SG(5) EG(5)];
G2f =   [Boron(6) HMS(6) AS(6) t300(6) Kevlar(6) SG(6) EG(6)];
v1f =   [Boron(7) HMS(7) AS(7) t300(7) Kevlar(7) SG(7) EG(7)];
v2f =   [Boron(8) HMS(8) AS(8) t300(8) Kevlar(8) SG(8) EG(8)];
SLpf =  [Boron(9) HMS(9) AS(9) t300(9) Kevlar(9) SG(9) EG(9)];
SLmf =  [Boron(10) HMS(10) AS(10) t300(10) Kevlar(10) SG(10) EG(10)];
ef1 =   [Boron(11) HMS(11) AS(11) t300(11) Kevlar(11) SG(11) EG(11)];

% Specify which material is being used in a variable
[F,OK] = listdlg('ListString',NameF);
[M,OK] = listdlg('ListString',NameM);

% Putting properties of chosen material into variables I can use
fd = fd(F);
pf = pf(F);
E1f = E1f(F);
E2f = E2f(F);
G1f = G1f(F);
G2f = G2f(F);
v1f = v1f(F);
v2f = v2f(M);
SLpf = SLpf(F)*1000;       % Have to convert from MPa to GPa
SLmf = SLmf(F)*1000;       % Have to convert from MPa to GPa
ef1 = ef1(F);

%%
% -----Matrix material-----
%   pm(g/in^3)   Em(GPa)   vm    S+(MPa)   S-(MPa)  SLT(MPa)   eTp
LM =        [1162556 2.21 .43 55 103 55 .081];         %LM
IMLS =      [1273276 3.45 .41 48 145 48 .014];       %IMLS
IMHS =      [1217916 3.45 .35 103.4 241.3 90 .02];   %IMHS
HM =        [1245596 5.17 .35 138 345 103 .02];        %HM
Polyimide = [1217916 3.45 .35 103 207 90 .02];  %Polyimide
PMR =       [1217916 3.24 .36 379 110 55 .02];        %PMR

% Organize the material properties
pm =        [LM(1) IMLS(1) IMHS(1) HM(1) Polyimide(1) PMR(1)];
Em =        [LM(2) IMLS(2) IMHS(2) HM(2) Polyimide(2) PMR(2)];
vm =        [LM(3) IMLS(3) IMHS(3) HM(3) Polyimide(3) PMR(3)];
SMP =     [LM(4) IMLS(4) IMHS(4) HM(4) Polyimide(4) PMR(4)];
SMM =    [LM(5) IMLS(5) IMHS(5) HM(5) Polyimide(5) PMR(5)];
SLT =       [LM(6) IMLS(6) IMHS(6) HM(6) Polyimide(6) PMR(6)];
eTp =       [LM(7) IMLS(7) IMHS(7) HM(7) Polyimide(7) PMR(7)];

% Putting properties of chosen material into variables I can use
pm = pm(M);
Em = Em(M);
vm = vm(M);
SMP = SMP(M);
SMM = SMM(M)*1000;
SLT = SLT(M)*1000;
eTp = eTp(M);
Gm = Em/(2*(1+vm));

% Find the volume fraction of the fibers and then the matrix
Vf = inputdlg('Volume fraction of the fiber (0<v_f<1):  ');
Vf = str2num(cell2mat(Vf));    % Converting the input string to a matrix
Vm = 1-Vf;                     % Notice that it's BIG 'V'm.

% Input the number of plys for the whole composite
n = inputdlg('Input the number of plys: ');
n = str2num(cell2mat(n));

%%

t_total = 0;    % Initializing total thickness for FOR loop
% For loop running through all the plys
T = zeros(3,3,n);
T2 = zeros(3,3,n);
for i = 1:n

% User input
prompt = {['Input the angle for ply ',num2str(i) ,' in degrees: '],['Input the thickness of ply ',num2str(i) ,' in mm: ']};
answer =  inputdlg(prompt);
thetatemp = str2num(cell2mat(answer(1)));
t(i) = str2num(cell2mat(answer(2)));
theta(i) = thetatemp*(pi/180);                 % Conversion to radians

t_total = t_total+t(i);               % Total thickness

% Total modulus of elasticity (Eqn 3.27)
E1 = E1f*Vf+Em*Vm;

% Major Poisson's Ratio of the fiber (Eqn 3.45)
v12 = v1f*Vf+vm*Vm;

% Effective transverse modulus (Eqn 3.54)
E2 = Em*((1-sqrt(Vf))+sqrt(Vf)/(1-sqrt(Vf)*(1-Em/E2f)));

% v21 (Eqn 2.25)
v21 = v12*(E2/E1);

% G1 (Same as E2 pretty much)
G1 = Gm*((1-sqrt(Vf))+sqrt(Vf)/(1-sqrt(Vf)*(1-Gm/G1f)));

% Equation to find the compliance matrix of this ply
S = [1/E1 -v12/E1 0;
    -v21/E2 1/E2 0;
    0 0 1/G1];

% Rather than doing Q = inv(S), I wrote out the formula for Q from page 68
% because MATLAB does their inverses in a slightly different way than we
% would.
Q = [E1/(1-v12*v21) v12*E2/(1-v12*v21) 0;
    v21*E1/(1-v12*v21) E2/(1-v12*v21) 0;
    0 0 G1];

Ttemp = [cos(theta(i))^2 sin(theta(i))^2 2*cos(theta(i))*sin(theta(i));
    sin(theta(i))^2 cos(theta(i))^2 -2*cos(theta(i))*sin(theta(i));
    -cos(theta(i))*sin(theta(i)) cos(theta(i))*sin(theta(i)) cos(theta(i))^2-sin(theta(i))^2];

T2temp = [cos(theta(i))^2 sin(theta(i))^2 cos(theta(i))*sin(theta(i));
    sin(theta(i))^2 cos(theta(i))^2 -cos(theta(i))*sin(theta(i));
    -2*cos(theta(i))*sin(theta(i)) 2*cos(theta(i))*sin(theta(i)) cos(theta(i))^2-sin(theta(i))^2];

% Getting the almost final Qbar and Sbar
Qbartemp = pinv(Ttemp)*Q*T2temp;
Sbartemp = pinv(T2temp)*S*Ttemp;

% Saving the T matrix into one total matrix
T(:,:,i) = Ttemp;
T2(:,:,i) = T2temp;

% % Setting up the 3 x 3 x i matrix for the input info
% SbarTotal = zeros(3,3,n);
% QbarTotal = zeros(3,3,n);

% Ending the for loop with filling in the 3 x 3 x i matrix
Sbar(:,:,i) = Sbartemp;
Qbar(:,:,i) = Qbartemp;

end

%%
% ---------------- Calculating z and ABBD matrix ----------------
% 't' still stands for thickness of each layer
t_total = sum(t(1:n));       % meters

% Setting up the z matrix
z = zeros(1,n+1);
z(1) = -t_total/2;      % half of the negative total thickness

% Setting up the A,B,D matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

% Finding the actual 'z's
for i = 2:n+1     % z position
    z(i) = z(1)+sum( t(1:i-1) );
    
     % Finding the A,B,D matrices with recently found 'z's
    A = A + Qbar(:,:,i-1).*(z(i)-z(i-1));
    B = B + Qbar(:,:,i-1).*.5*(z(i)^2-z(i-1)^2);
    D = D + Qbar(:,:,i-1).*1/3*(z(i)^3-z(i-1)^3);
end

% Setting up the ABBD matrix and taking inverse
ABBD = [A B;
    B D];
invABBD = pinv(ABBD);

%%
% Input the Nx, Ny, Nxy, Mx, My, and Mxy numbers
prompt = {'Input Nx in GPa-mm: ','Input Ny in GPa-mm: ','Input Nxy in GPa-mm: ','Input the Mx in GPa-mm^2: ','Input the My in GPa-mm^2: ','Input the Mxy in GPa-mm^2: '};
answer = inputdlg(prompt);
NM = str2num(cell2mat(answer'));      % str2double also works

%%
% Finding the Eps0 and k(appa) arrays
Deformations = invABBD*NM;
Eps0 = Deformations(1:3);     % Epsilon0
k = Deformations(4:6);        % Kappa (Curvatures)

% Finding the 3xn matrix for each part of x, y, and xy with epsilon and sigma.
epsbarT = zeros(3,n);
epsbarB = zeros(3,n);
for i = 1:3
    for j = 1:n
        epsbarT(i,j) = Eps0(i) +z(j)*k(i);
        epsbarB(i,j) = Eps0(i) +z(j+1)*k(i);
    end
end

% Finding the sigma in x, y, and xy directions for each layer.
for i = 1:n
    sigmabarT(:,i) = Qbar(:,:,i)*epsbarT(:,i);
    sigmabarB(:,i) = Qbar(:,:,i)*epsbarB(:,i);
end

%% 
% Find the stresses and strains on the top and bottom of each ply in 1-2
% frame
for i = 1:n
    sigma12T = T(i)*sigmabarT(:,i);
    sigma12B = T(i)*sigmabarB(:,i);
    
    eps12T = T(i)*epsbarT(:,i);
    eps12B = T(i)*epsbarB(:,i);
end

% Apply relevant failure model
% Find SLp
Smf1p = Em*ef1;     % Find the stress that the matrix is at when the fiber fails
SLP = SLpf*Vf + Smf1p*Vm;     % Eqn -

% Find STp
DoverS = sqrt(4*Vf/pi);       % In slide in class notes before Eqn 4.34
F = 1/(DoverS*(Em/E2f-1)+1);   % Eqn 4.38 in class notes
eTp = eTp/F;                  % Eqn 4.34 in class notes
STP = E2*SMP/(Em*F);         % Eqn 4.36 in book

% Find SLm
SLM = E1*eTp/v12;     % Eqn 4.33

% Find STm (matrix Sminus)
STM = SMM;

% Find SLT (matrix SLT)
SLT;

%% 
% Have the user choose a failure method
failure = {'Tsai-Hill','Max Stress','Max Strain'};
[X,OK] = listdlg('ListString',failure);

% Output list of layers that fail
% Pull out the Sigmas and Taus for the top and bottom of each layer
    sigma1T = sigmabarT(1,:);
    sigma2T = sigmabarT(2,:);
    tau12T = sigmabarT(3,:);
    sigma1B = sigmabarB(1,:);
    sigma2B = sigmabarB(2,:);
    tau12B = sigmabarB(3,:);
%%%%%%%%%%%%%%%%%%%%%% Tsai-Hill Method %%%%%%%%%%%%%%%%%%%%%%
if X == 1;
    % Run through the top of each layer
    for i = 1:n   
        % If sigma1 is positive, use SLP. If negative, use SLM
        if sigma1T(1,i)>0
            SL(i) = SLP;
        else
            SL(i) = SLM;
        end

        % If sigma2 is positive, use STP. If negative, use STM
        if sigma1T(1,i)>0
            ST(i) = SLP;
        else
            ST(i) = SLM;
        end

        % Calculate the Tsai-Hill Value for each layer
        TsaiB(i) = sigma1T(1,i)^2/SL(i)^2-sigma1T(1,i)*sigma2T(1,i)/SL(i)^2+sigma2T(1,i)^2/ST(i)^2+tau12T(1,i)^2/SLT^2;
    end
    % See which top layers fail
    for i = 1:n
        if TsaiB(i) > 1
            output = ['[The top of layer ',num2str(i),' failed.]'];
            disp(output)
        else
            output = ['[The top of layer ',num2str(i),' did not fail.]'];
            disp(output)
        end
    end
    
    % Run through the bottom of each layer
    for i = 1:n   
        % If sigma1 is positive, use SLP. If negative, use SLM
        if sigma1B(1,i)>0
            SL(i) = SLP;
        else
            SL(i) = SLM;
        end

        % If sigma2 is positive, use STP. If negative, use STM
        if sigma1B(1,i)>0
            ST(i) = SLP;
        else
            ST(i) = SLM;
        end

        % Calculate the Tsai-Hill Value for each layer
        TsaiB(i) = sigma1B(1,i)^2/SL(i)^2-sigma1B(1,i)*sigma2B(1,i)/SL(i)^2+sigma2B(1,i)^2/ST(i)^2+tau12B(1,i)^2/SLT^2;
    end
    % See which bottom layers fail
    for i = 1:n
        if TsaiB(i) > 1
            output = ['[The bottom of layer ',num2str(i),' failed.]'];
            disp(output)
        else
            output = ['[The bottom of layer ',num2str(i),' did not fail.]'];
            disp(output)
        end
    end

%%%%%%%%%%%%%%%%%%%%%% Max Stress Method %%%%%%%%%%%%%%%%%%%%%%
elseif X == 2;
    % Run through the top of each layer
    for i = 1:n   
        % If sigma1 is positive, use SLP. If negative, use SLM
        if sigma1T(1,i)>0
            SL(i) = SLP;
        else
            SL(i) = SLM;
        end

        % If sigma2 is positive, use STP. If negative, use STM
        if sigma1T(1,i)>0
            ST(i) = SLP;
        else
            ST(i) = SLM;
        end
        
        % See which top layer fails
        if abs(sigma1T(i)) > SL(i) || abs(sigma2T(i)) > ST(i) || abs(tau12T(i)) > SLT
            output = ['[The top of layer ',num2str(i),' failed.]'];
            disp(output)
        else
            output = ['[The top of layer ',num2str(i),' did not fail.]'];
            disp(output)
        end
        
        % See which bottom layer fails
        if abs(sigma1B(i)) > SL(i) || abs(sigma2B(i)) > ST(i) || abs(tau12B(i)) > SLT
            output = ['[The bottom of layer ',num2str(i),' failed.]'];
            disp(output)
        else
            output = ['[The bottom of layer ',num2str(i),' did not fail.]'];
            disp(output)
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%% Max Strain Method %%%%%%%%%%%%%%%%%%%%%%
else X == 3;
    for i = 1:n
    % Define slopes in Max Strain Parallelogram
    Slope1 = 1/v12;
    Slope2 = v21;
    
    % Y-value of sigma1 
    ySTP(i) = Slope1*sigma1T(i) + STP;
    ySTM(i) = Slope1*sigma1T(i) + STM;
    % Use point slope form [y-y1 = m(x-x1)] for the vertical plot lines
    ySLP(i) = Slope2*sigma1T(i)-SLP;     % y1 = 0, x1 = x-intercept
    ySLM(i) = Slope2*sigma1T(i)-SLM;     % y1 = 0, x1 = x-intercept
    
        % See which top layer fails
        if abs(sigma2T(i)) > ySTP(i) || abs(sigma2T(i)) > abs(ySTM(i)) || abs(tau12T(i)) > SLT
            output = ['[The top of layer ',num2str(i),' failed.]'];
            disp(output)
        else
            output = ['[The top of layer ',num2str(i),' did not fail.]'];
            disp(output)
        end
        
        % See which bottom layer fails
        if abs(sigma2B(i)) > ySTP(i) || abs(sigma2B(i)) > abs(ySTM(i)) || abs(tau12B(i)) > SLT
            output = ['[The top of layer ',num2str(i),' failed.]'];
            disp(output)
        else
            output = ['[The top of layer ',num2str(i),' did not fail.]'];
            disp(output)
        end
    end
end

% Done.