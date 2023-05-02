%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of T3                          *
%*                                                                 *
%* - Call procedures:                                              *
%*     T3Stiff.m - InitT3()                                    *
%*     ReadT3.m - ReadT3()                                     *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     Yourself                                                    *
%*                                                                 *
%* *****************************************************************

function T3Stiff()

% Init variables of the element
InitT3();

% Read Material and Elements
ReadT3();

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
function InitT3()
global sdata;
sdata.NNODE = 3;
sdata.NDOF = 2;

end

% Assemble structure stiffness matrix
function Assemble()
% To be completed

end

