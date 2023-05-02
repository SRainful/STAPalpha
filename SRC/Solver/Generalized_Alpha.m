%* ************************************************************************
%* - Function of STAPMAT for Time Integral in Generalized-Alpha Method    *
%*                                                                        *
%* - Purpose:                                                             *
%*                    *
%*                                                                        *
%* - Call procedures:                                                     *
%*                                                                        *
%* - Called by :                                                          *
%*     Solve.m                                                          *
%*                                                                        *
%* - Programmed by:                                                       *
%*     TianYu Zhao                                                        *              
%*                                                                        *
%* ************************************************************************
function Generalized_Alpha(L)

global cdata;
global sdata;

%Get Mass Matrix
GetMass();

%Get Damp Matrix
GetDamp(L);

%Get Initial position
InitialPos();

%Time Integral
TimeIntegral(L);

end

% ----------------------- Functions -----------------------------------
function GetMass()

global cdata;
global sdata;

if ~(sdata.MASS == 0) 
% 如果质量阵已经在之前工况被计算过，不再计算
    return
else
NUME = cdata.NPAR(2);           %单元个数
AREA = sdata.AREA;              %截面面积
RHO = sdata.RHO;                %密度
XYZ = sdata.XYZ;                %节点位置
MATP = sdata.MATP;              %单元对应的材料属性编号

NEQ = sdata.NEQ;
LM = sdata.LM;
NDOF = sdata.NDOF;

%Initial concentrate mass matrix
MASS = zeros(NEQ , 1);

for N = 1 : NUME
%   Calculate the length of element
    DX = XYZ(1,N) - XYZ(4,N);
    DY = XYZ(2,N) - XYZ(5,N);
    DZ = XYZ(3,N) - XYZ(6,N);
    DL = sqrt(DX^2 + DY^2 + DZ^2);

%   Calculate the mass of element
    MTYPE = MATP(N);
    MRHO = RHO(MTYPE);
    MAREA = AREA(MTYPE);
    W = MRHO * MAREA * DL;

%   Assemble concentrate mass matrix
    for I = 1:NDOF
        F = LM(I,N);
        if (F > 0)
            MASS(F) = MASS(F) + W/2;
        end
    end
    for J = 1+NDOF : 3+NDOF
        F = LM(J,N);
        if (F > 0)
            MASS(F) = MASS(F) + W/2;
        end
    end
end

sdata.MASS = sparse(diag(MASS));

end

end


function GetDamp(L)

global sdata;
global cdata;
ALPHA = cdata.ALPHA;
BETA = cdata.BETA;

% Rayleigh Damp
DAMP = ALPHA(L) * sdata.MASS + BETA(L) * sdata.STIFFOrigin;
sdata.DAMP = DAMP;
C = full(DAMP);

end

function InitialPos()

global sdata; 
global cdata;

ID = sdata.ID;
NUMNP = cdata.NUMNP;
POS = zeros(sdata.NEQ,1);
X = sdata.X;
Y = sdata.Y;
Z = sdata.Z;

for I = 1:3
    for J = 1:NUMNP
        F= ID(I,J);
        if(F > 0)
            switch I
                case 1
                    POS(F) = X(J);
                case 2
                    POS(F) = Y(J);
                case 3
                    POS(F) = Z(J);
            end
        end
    end
end

sdata.POS = POS;

end


function TimeIntegral(L)

global cdata;
global sdata;

% Initial STIFF,MASS,DAMP matrix
K = sdata.STIFFOrigin;
M = sdata.MASS;
C = sdata.DAMP;

% Load conditions
Q = sdata.Q;

% Initial conditions
V = sdata.V;
POS = sdata.POS;
X = zeros(sdata.NEQ,1);
A = M\(Q-K*X-C*V);


% Integral parameters
T = cdata.MTIME(L);
Dt = cdata.MDELTAT(L);
STEPT = 0:Dt:T;

% Parameters for Generalized Alpha method calculating 
R_inf = 1;
Af = R_inf/(R_inf+1);
Am = (2*R_inf-1)/(R_inf+1);
beta = (1-Am+Af)^2/4;
gamma = 0.5-Am+Af;
    
ck = 1-Af;
c0 = (1-Am)/(beta*Dt^2);
c1 = (ck*gamma)/(beta*Dt);
c2 = Dt*c0;
c3 = c2*Dt/2-1;
c4 = ck*gamma/beta-1;
c5 = ck*(gamma/(2*beta)-1)*Dt;

% Integral process
cdata.IDYNA = fopen('.\Data\DYNA.OUT', 'w');
IDYNA = cdata.IDYNA;
ID = sdata.ID;
XI = sdata.X;
YI = sdata.Y;
ZI = sdata.Z;

fprintf(IDYNA, 'T I M E   A N D   N O D A L   P O S I T I O N   D A T A \n\n');
fprintf(IDYNA, '       TIME       NODE NUMBER             POSITION\n');
fprintf(IDYNA, '                                   X          Y          Z\n');

for I = 1:length(STEPT)

    D = zeros(3, 1, 'double');
    for NOD = 1 : cdata.NUMNP
        D(:) = 0;
        if (ID(1, NOD) ~= 0) 
            D(1) = POS(ID(1, NOD)) + X(ID(1, NOD));
        else
            D(1) = XI(NOD);
        end

        if (ID(2, NOD) ~= 0) 
            D(2) = POS(ID(2, NOD)) + X(ID(2, NOD));
        else
            D(2) = YI(NOD);
        end

        if (ID(3, NOD) ~= 0) 
            D(3) = POS(ID(3, NOD)) + X(ID(3, NOD));
        else
            D(3) = ZI(NOD);
        end
        fprintf(IDYNA, '       %.3f      %6d      %10.3f%10.3f%10.3f\n' ,...
                STEPT(I),NOD,D(1),D(2),D(3));
    end

    Kb = ck*K + c0*M + c1*C;
    Qb = Q - Af*K*A + M*(c0*X + c2*V + c3*A) + C*(c1*X + c4*V + c5*A);
    % t + Δt 
    Xw = Kb \ Qb;
    Vw = gamma/(beta*Dt)*(Xw - X) + (1-gamma/beta)*V + (1-gamma/(2*beta))*Dt*A;
    Aw = 1/(beta*Dt^2)*(Xw - X) - 1/(beta*Dt)*V - (1/(2*beta)-1)*A;
    
    
    X = Xw;
    V = Vw;
    A = Aw;

end


end



