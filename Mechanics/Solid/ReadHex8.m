%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of Solid3d                     *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadSolid.m - ReadMaterial()                                *
%*     ReadSolid.m - ReadElements()                                *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Hex8Stiff.m                                              *
%*                                                                 *
%* - Programmed by:                                                *
%*     Yourself                                                    *
%*                                                                 *
%* *****************************************************************

function ReadHex8()

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
% To be completed

end

% Read elements information
function ReadElements()
% To be completed

end