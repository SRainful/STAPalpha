%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of Beam                        *
%*                                                                 *
%* - Call procedures:                                              *
%*     BeamStiff.m - InitBeam()                                    *
%*     ReadBeam.m - ReadBeam()                                     *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     Yourself                                                    *
%*                                                                 *
%* *****************************************************************

function BeamStiff()

% Init variables of the element
InitBeam();

% Read Material and Elements
ReadBeam();

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres();

% Data check Or Solve
global cdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble();

end

% ----------------------- Functions -----------------------------------

% Init parameters of beam element
function InitBeam()
global sdata;
sdata.NNODE = 2;
sdata.NDOF = 6;  % For Euler-Bernoulli Beam

end

% Assemble structure stiffness matrix
function Assemble()
% To be completed
% Get global data
global sdata;
global cdata;

sdata.STIFF = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
E = sdata.E; AREA = sdata.AREA; Iy = sdata.Iy; Iz = sdata.Iz; Jx = sdata.Jx; NU = sdata.NU; LM = sdata.LM;
G = E./(2*(1+NU));



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
    
%   SRC/Mechanics/ADDBAN.m
    ADDBAN(S, LM(:, N));
    
end

% The third time stamp
cdata.TIM(3, :) = clock;

end


