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

% Init parameters of truss element
function InitBeam()
global sdata;
sdata.NNODE = 2;
sdata.NDOF = 6;  % For Euler-Bernoulli Beam

end

% Assemble structure stiffness matrix
function Assemble()
% To be completed

end


