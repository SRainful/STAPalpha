%* ************************************************************************
%* - Function of STAPMAT for Time Integral in Generalized-Alpha Method    *
%*                                                                        *
%* - Purpose:                                                             *
%*                    *
%*                                                                        *
%* - Call procedures:                                                     *
%*                                                                        *
%* - Called by :                                                          *
%*     stapmat.m                                                          *
%*                                                                        *
%* - Programmed by:                                                       *
%*     TianYu Zhao                                                        *              
%*                                                                        *
%* ************************************************************************
function Generalized_Alpha();

global cdata;
global sdata;

%Get Mass Matrix
GetMass();

%Get Damp Matrix
GetDamp(L);

%Time Integral
TimeIntegral(L);

end

% ----------------------- Functions -----------------------------------
function GetMass();

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
NODEOFELE = sdata.NODEOFELE;    %单元对应的节点编号
MATP = sdata.MATP;              %单元对应的材料属性编号

NUMNP = cdata.NUMNP;
NDOF = sdata.NDOF;

%Initial concentrate mass matrix
MASS = zeros(NUMNP * NDOF , 1);

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
    RANI = NDOF * (NODEOFELE(1,N) - 1) + 1 : 1 : NDOF * (NODEOFELE(1,N) - 1) + 3;
    RANJ = NDOF * (NODEOFELE(2,N) - 1) + 1 : 1 : NDOF * (NODEOFELE(2,N) - 1) + 3;
    MASS(RANI) = MASS(RANI) + 0.5 * W * ones(3,1);    
    MASS(RANJ) = MASS(RANJ) + 0.5 * W * ones(3,1);     

end

sdata.MASS = MASS;

end

end


function GetDamp(L)

global sdata;

% Rayleigh Damp
DAMP = ALPHA(L) * sdata.MASS + BETA(L) * sdata.STIFF;
sdata.DAMP = DAMP;

end


function TimeIntegral(L)

global cdata;
global sdata;




end



