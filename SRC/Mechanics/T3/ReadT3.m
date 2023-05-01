%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of T3                          *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadT3.m - ReadMaterial()                                *
%*     ReadT3.m - ReadElements()                                *
%*                                                                 *
%* - Called by :                                                   *
%*     ./T3Stiff.m                                              *
%*                                                                 *
%* - Programmed by:                                                *
%*     Yourself                                                    *
%*                                                                 *
%* *****************************************************************

function ReadT3()

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