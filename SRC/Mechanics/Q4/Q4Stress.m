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
%*     Zhao TY                                                    *
%*                                                                 *
%* *****************************************************************

function Q4Stress(NUM, NG)%载荷工况编号/节点自由度编号
% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
NUME = sdata.NUME; 
MATP = sdata.MATP; 
XYZ  = sdata.XYZ;
E    = sdata.E; 
NU   = sdata.NU; 
LM   = sdata.LM;
CON  = cdata.NPAR(1);
U    = sdata.DIS(:, NUM);%有自由度的节点上位移

fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  O N  G A U S S  P O I N T 1  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT             EPSILON-X           EPSILON-Y           EPSILON-XY            SIGMA-X            SIGMA-Y             SIGMA-XY \n' ...
    '       NUMBER\n'], NG);
% For every element
    % Calculate physical matrix
    % Compute the strain and stress matrix on Gauss points
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
    %Get displacement of nodes
    di = zeros(8,1);
    step = 1;

    for k = 1:8
        if LM(k,N)>0
            di(k) = U(step);
            step = step+1;
        else
            di(k) = double(0);
        end
    end

    % Compute the strain and stress matrix on Gauss points
    % Set Guass integral point
    GP = [-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3);...
          1/sqrt(3),-1/sqrt(3),-1/sqrt(3),1/sqrt(3)];
    % Compute Jacobi matrix's det
    DetJ = zeros(1,4,'double'); 
    XY = [XYZ(1,N),XYZ(2,N);...
          XYZ(3,N),XYZ(4,N);...
          XYZ(5,N),XYZ(6,N);...
          XYZ(7,N),XYZ(8,N);];
   
    STRESS = zeros(3,4);
    STRAIN = zeros(3,4);
    for i = 1:4
        neta = GP(2,i);
        epsi = GP(1,i);       
        Jaci = 1/4*[neta-1,1-neta,1+neta,-neta-1;...
                    epsi-1,-epsi-1,1+epsi,1-epsi;]*XY;
        DetJ(i) = det(Jaci);    
    % Compute B      
        Jxy = inv(Jaci)*1/4*[neta-1,1-neta,1+neta,-neta-1;...
                             epsi-1,-epsi-1,1+epsi,1-epsi;];
        Bi = zeros(3,8);
        Bi(1,[1,3,5,7]) = Jxy(1,:);
        Bi(2,[2,4,6,8]) = Jxy(2,:);
        Bi(3,[1,3,5,7]) = Jxy(2,:);
        Bi(3,[2,4,6,8]) = Jxy(1,:);
        STRAIN(:,i) = Bi*di;%应变
        STRESS(:,i) = D*STRAIN(:,i);%应力
    end
  
    fprintf(IOUT, ' %10d             %13.6e       %13.6e       %13.6e       %13.6e       %13.6e       %13.6e\n', N,  STRAIN(3,1),STRAIN(3,2),STRAIN(3,3), STRESS(3,1),STRESS(3,2),STRESS(3,3));
end    
% NOTES: The stress messages are on Gauss point
% You'd better map them on the nodes. 
% It's better if you can output the results in the TECPLOT form for a wonderful visualization.

end
