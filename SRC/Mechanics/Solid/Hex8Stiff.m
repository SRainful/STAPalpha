%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of Hex8                        *
%*                                                                 *
%* - Call procedures:                                              *
%*     Hex8Stiff.m - InitHex8()                                    *
%*     ReadHex8.m - ReadHex8()                                     *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     Yourself                                                    *
%*                                                                 *
%* *****************************************************************

function Hex8Stiff()

% Init variables of the element
InitHex8();

% Read Material and Elements
ReadHex8();

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
function InitHex8()
global sdata;
sdata.NNODE = 8;
sdata.NDOF = 3;

end

% Assemble structure stiffness matrix
function Assemble()
% To be completed

end


