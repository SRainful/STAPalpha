%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of Q4                          *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadQ4.m - ReadMaterial()                                   *
%*     ReadQ4.m - ReadElements()                                   *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Q4Stiff.m                                                 *
%*                                                                 *
%* - Programmed by:                                                *
%*     Zhao TY                                                    *
%*                                                                 *
%* *****************************************************************

function ReadQ4()

% Read Material information
ReadMaterial()

% Read Element information
ReadElements()

% the second time stamp
global cdata;
cdata.TIM(2,:) = clock;

end

% ----------------------- Functions -----------------------------------
% Read Material information
function ReadMaterial()
% Get global data
global cdata;
global sdata;
% Get file pointers
IIN = cdata.IIN;
IOUT = cdata.IOUT;
% Write Echo message here
if (cdata.NPAR(3) == 0) cdata.NPAR(3) = 1; end
fprintf(IOUT, '\n M A T E R I A L   D E F I N I T I O N\n');
fprintf(IOUT, '\n NUMBER OF DIFFERENT SETS OF MATERIAL\n');
fprintf(IOUT, ' AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . = %10d\n', cdata.NPAR(3));
fprintf(IOUT, '  SET       YOUNG''S     CROSS-SECTIONAL\n');
fprintf(IOUT, ' NUMBER     MODULUS          NU\n');
fprintf(IOUT, '               E              NU\n');
% Read material message here
sdata.NUME = cdata.NPAR(2);
sdata.NUMMAT = cdata.NPAR(3);
NUMMAT = cdata.NPAR(3);%材料种类数
sdata.E = zeros(NUMMAT, 1);%刚度
sdata.NU = zeros(NUMMAT, 1);%泊松比
for I = 1:cdata.NPAR(3)
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E(N) = tmp(2);
    sdata.NU(N) = tmp(3);
    fprintf(IOUT, '%5d    %12.5e  %14.6e\n', N, tmp(2), tmp(3));
end

end

% Read elements information
function ReadElements()
% Get global data
global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
IOUT = cdata.IOUT;

% Write Echo message here

fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n');
fprintf(IOUT, '\n      ELEMENT          NODE          NODE          NODE          NODE       MATERIAL\n');
fprintf(IOUT, '      NUMBER-N           I             J             K             M        SET NUMBER\n');
% Get Position data here
NUME = cdata.NPAR(2);
sdata.XYZ = zeros(2*4, NUME);
sdata.MATP = zeros(NUME, 1);     
sdata.LM = zeros(2*4, NUME );    % Connectivity Matrix
sdata.MHT = zeros(sdata.NEQ, 1);
X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;
XYZ = sdata.XYZ; MATP = sdata.MATP; LM = sdata.LM;

for N = 1:NUME
    % The sequence of nodal ID in anti-clockwise
    tmp = str2num(fgetl(IIN));
    I = round(tmp(2));
    J = round(tmp(3));
    K = round(tmp(4)); 
    M = round(tmp(5));    
    MTYPE = round(tmp(6));
    
%   Save element information
    XYZ(1, N) = X(I);
    XYZ(2, N) = Y(I);
    XYZ(3, N) = X(J);
    XYZ(4, N) = Y(J);
    XYZ(5, N) = X(K);
    XYZ(6, N) = Y(K);
    XYZ(7, N) = X(M);
    XYZ(8, N) = Y(M);  
    MATP(N) = MTYPE;
    
    fprintf(IOUT, '%10d      %10d    %10d    %10d    %10d       %5d\n', N, I, J,K,M, MTYPE);

%   Compute connectivity matrix
    LM(1, N) = ID(1, I);
    LM(3, N) = ID(1, J);
    LM(5, N) = ID(1, K);
    LM(7, N) = ID(1, M);    
    LM(2, N) = ID(2, I);
    LM(4, N) = ID(2, J);
    LM(6, N) = ID(2, K);
    LM(8, N) = ID(2, M);    

%   Updata column heights and bandwidth
    ColHt(LM(:, N))
end
sdata.XYZ = XYZ; sdata.MATP = MATP; sdata.LM = LM;
% Clear the memory of X, Y, Z
sdata.X = double(0);
sdata.Y = double(0);
sdata.Z = double(0);
end