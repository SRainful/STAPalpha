%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of Q4                          *
%*                                                                 *
%* - Call procedures:                                              *
%*     Q4Stiff.m - InitQ4()                                        *
%*     ReadQ4.m - ReadQ4()                                         *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     Zhao TY                                                    *
%*                                                                 *
%* *****************************************************************

function Q4Stiff()

% Init variables of the element
InitQ4();

% Read Material and Elements
ReadQ4();

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

% Init parameters of truss element
function InitQ4()
    global sdata;
    sdata.NNODE = 4;
    sdata.NDOF = 2;
end

% Assemble structure stiffness matrix
function Assemble()
% Get global data
global sdata; 
global cdata; 

sdata.STIFF = zeros(sdata.NWK, 1, 'double');%总刚

NUME = sdata.NUME; 
MATP = sdata.MATP; 
XYZ = sdata.XYZ;
E = sdata.E;
NU = sdata.NU; 
LM = sdata.LM;
CON=cdata.NPAR(1);
% For every elements
    % Calculate physical matrix
    % Set Guass integral point
    % Compute N,B,K,Jacobi matrix at Guass point
    % Set Guass integral point
 %四点高斯积分
GP = [-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3);...
      1/sqrt(3),-1/sqrt(3),-1/sqrt(3),1/sqrt(3)];
W = ones(1,4);
for N = 1:NUME
    MTYPE = MATP(N);
    % Calculate physical matrix
    if CON == 2
    D = E(MTYPE) /(1- NU(MTYPE)^2)*...
        [1,NU(MTYPE),0;...
         NU(MTYPE),1,0;...
         0,0,(1-NU(MTYPE))/2];
    else
    D = E(MTYPE) /((1+ NU(MTYPE))*(1-2*NU(MTYPE)))*...
        [1,NU(MTYPE),0;...
         NU(MTYPE),1,0;...
         0,0,(1-2*NU(MTYPE))/2];    
    end

    % Compute Jacobi matrix's det
    K = zeros(8,8,'double');%单刚

    XY = [XYZ(1,N),XYZ(2,N);...
          XYZ(3,N),XYZ(4,N);...
          XYZ(5,N),XYZ(6,N);...
          XYZ(7,N),XYZ(8,N);];
   
    for i = 1:4
        neta = GP(2,i);
        epsi = GP(1,i);       
        Jaci = 1/4*[neta-1,1-neta,1+neta,-neta-1;...
                    epsi-1,-epsi-1,1+epsi,1-epsi;]*XY;
        DetJi = det(Jaci);    
    % Compute B      
        Jxy = inv(Jaci)*1/4*[neta-1,1-neta,1+neta,-neta-1;...
                             epsi-1,-epsi-1,1+epsi,1-epsi;];
        Bi = zeros(3,8);
        Bi(1,[1,3,5,7]) = Jxy(1,:);
        Bi(2,[2,4,6,8]) = Jxy(2,:);
        Bi(3,[1,3,5,7]) = Jxy(2,:);
        Bi(3,[2,4,6,8]) = Jxy(1,:);
    % Compute K
        K  = K+Bi'*D*Bi*DetJi;
    end   
    % Save K in LM with function ADDBAN()
    ADDBAN(K, LM(:, N));
end
cdata.TIM(3, :) = clock;
end


