%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* - Programmed by:                                                *
%*     Yourself                                                    *
%*                                                                 *
%* *****************************************************************

function BeamStress(NUM, NG)
% To be completed
% Get global data


global cdata;
global sdata;
global fname;

f_name = strcat('.\Data\', fname);
cdata.IIN = fopen(f_name, 'r');
IIN = cdata.IIN;
IOUT = cdata.IOUT;
NUMNP=cdata.NUMNP;NLOAD=cdata.NLOAD;NODEOFELE = sdata.NODEOFELE;NUME = sdata.NUME; MATP = sdata.MATP;NUMMAT=sdata.NUMMAT;XYZ = sdata.XYZ;
E = sdata.E; AREA = sdata.AREA; Iy = sdata.Iy; Iz = sdata.Iz; Jx = sdata.Jx; NU = sdata.NU; LM = sdata.LM;
G = E./(2*(1+NU));
U = sdata.DIS(:, NUM);

fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT             FORCE_X1        FORCE_Y1        FORCE_Z1        MOMENT_X1        MOMENT_Y1        MOMENT_Z1        FORCE_X2        FORCE_Y2        FORCE_Z2        MOMENT_X2        MOMENT_Y2        MOMENT_Z2\n' ...
    '       NUMBER\n'], NG);

GaussianCollection = zeros(3,NUME);
StressCollection = zeros(6,NUME);
% For every element
for N = 1:NUME
    MTYPE = MATP(N);
   
%   compute the length of beam element
    DX = XYZ(1, N) - XYZ(4, N);
    DY = XYZ(2, N) - XYZ(5, N);
    DZ = XYZ(3, N) - XYZ(6, N);
    XL2 = DX*DX + DY*DY + DZ*DZ;
    XL = sqrt(XL2);
    
%   坐标转换矩阵
    T = zeros(12, 12);
    CXX = (XYZ(4, N)-XYZ(1, N))/XL;
    CXY = (XYZ(5, N)-XYZ(2, N))/XL;
    CXZ = (XYZ(6, N)-XYZ(3, N))/XL;
    if (abs(CXX-0.0) < 10E-5) && (abs(CXY-0.0) < 10E-5)
        CYX = 1.0;
        CYY = 0.0;
        CYZ = 0.0;
        CZX = 0.0;
        CZY = 1.0;
        CZZ = 0.0;
    
    else
        CYX = -CXX*CXZ/sqrt(CXX*CXX+CXY*CXY);
        CYY = -CXY*CXZ/sqrt(CXX*CXX+CXY*CXY);
        CYZ = sqrt(CXX*CXX+CXY*CXY);
        CZX = CXY/sqrt(CXX*CXX+CXY*CXY);
        CZY = -CXX/sqrt(CXX*CXX+CXY*CXY);
        CZZ = 0.0;
    end
    T(1,1) = CXX;
    T(1,2) = CXY;
    T(1,3) = CXZ;
    T(2,1) = CYX;
    T(2,2) = CYY;
    T(2,3) = CYZ;
    T(3,1) = CZX;
    T(3,2) = CZY;
    T(3,3) = CZZ;
    for i = 1:3
        for j = 1:3
            T(i+3,j+3) = T(i,j);
            T(i+6,j+6) = T(i,j);
            T(i+9,j+9) = T(i,j);
        end
    end
%   计算刚度阵
    S = zeros(12, 12);
    S(1,1) = E(MTYPE)*AREA(MTYPE)/XL;
    S(2,2) = 12.0*E(MTYPE)*Iz(MTYPE)/(XL*XL*XL);
    S(3,3) = 12.0*E(MTYPE)*Iy(MTYPE)/(XL*XL*XL);
    S(4,4) = G(MTYPE)*Jx(MTYPE)/XL;
    S(5,5) = 4.0*E(MTYPE)*Iy(MTYPE)/XL;
    S(6,6) = 4.0*E(MTYPE)*Iz(MTYPE)/XL;
    S(7,7) = E(MTYPE)*AREA(MTYPE)/XL;
    S(8,8) = 12.0*E(MTYPE)*Iz(MTYPE)/(XL*XL*XL);
    S(9,9) = 12.0*E(MTYPE)*Iy(MTYPE)/(XL*XL*XL);
    S(10,10) = G(MTYPE)*Jx(MTYPE)/XL;
    S(11,11) = 4.0*E(MTYPE)*Iy(MTYPE)/XL;
    S(12,12) = 4.0*E(MTYPE)*Iz(MTYPE)/XL;
    S(5,3) = -6.0*E(MTYPE)*Iy(MTYPE)/(XL*XL);
    S(6,2) = 6.0*E(MTYPE)*Iz(MTYPE)/(XL*XL);
    S(7,1) = -E(MTYPE)*AREA(MTYPE)/XL;
    S(8,2) = -12.0*E(MTYPE)*Iz(MTYPE)/(XL*XL*XL);
    S(8,6) = -6.0*E(MTYPE)*Iz(MTYPE)/(XL*XL);
    S(9,3) = -12.0*E(MTYPE)*Iy(MTYPE)/(XL*XL*XL);
    S(9,5) = 6.0*E(MTYPE)*Iy(MTYPE)/(XL*XL);
    S(10,4) = -G(MTYPE)*Jx(MTYPE)/XL;
    S(11,3) = -6.0*E(MTYPE)*Iy(MTYPE)/(XL*XL);
    S(11,5) = 2.0*E(MTYPE)*Iy(MTYPE)/XL;
    S(11,9) = 6.0*E(MTYPE)*Iy(MTYPE)/(XL*XL);
    S(12,2) = 6.0*E(MTYPE)*Iz(MTYPE)/(XL*XL);
    S(12,6) = 2.0*E(MTYPE)*Iz(MTYPE)/XL;
    S(12,8) = -6.0*E(MTYPE)*Iz(MTYPE)/(XL*XL);
    
    for i = 1:11
        for j = i+1:12
            S(i,j) = S(j,i);
        end
    end
    
    S = T'*S*T;
    
%   计算每个单元的位移
    UELE = zeros(1,12);
    for i = 1:12
        if LM(i, N) > 0
            UELE(i) = U(LM(i, N));
        end
    end
   
    FORCE = zeros(1,12);
%   计算单元内力和力矩
    for i = 1:12
        Force_initial = 0;
        for j = 1:12
            Force_initial = Force_initial+UELE(j)*S(i,j); 
        end
        FORCE(i) = Force_initial;
    end
    
    fprintf(IOUT, ' %10d           %13.6e    %13.6e    %13.6e    %13.6e    %13.6e    %13.6e    %13.6e    %13.6e    %13.6e    %13.6e    %13.6e    %13.6e\n', N, FORCE(1), FORCE(2), FORCE(3), FORCE(4), FORCE(5), FORCE(6), FORCE(7), FORCE(8), FORCE(9), FORCE(10), FORCE(11), FORCE(12));
%     GaussianCollection(:,N) = 0.5*(XYZ(4:6,N)+XYZ(1:3,N));
    %StressCollection([1,2,3,5,6,4],N) = 0.5*(FORCE(7:12)-FORCE(1:6))/AREA(MATP(N));
    StressCollection([1,2,3,5,6,4],N) = 0.5*(FORCE(7:12)-FORCE(1:6));
end
% It's better if you can output the results in the TECPLOT form for a wonderful visualization.
NODEOFELE = NODEOFELE';
PostProcessor(NODEOFELE, 1, StressCollection);%（NODEOFELE可能要存）
end