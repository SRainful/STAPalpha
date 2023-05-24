% ******************************************************************
%*                                                                 *
%*                         S T A P M A T                           *
%*                                                                 *
%*      An In-CORE SOLUTION STATIC ANALYSIS PROGRAM IN MATLAB      *
%*      Adapted from STAP90 (FORTRAN 90) for teaching purpose      *
%*                                                                 *                              *
%*  - Programmed by:                                               *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%*     Modified version by Wang Shuai    2021.01.06                *
%*                                                                 *
%* *****************************************************************

% Set paths of functions
AddPath();

% Define Global Variables
global cdata;
global sdata;
global fname;

cdata = ControlData;
sdata = SolutionData;

% Read InPut file
fname = 'Beam3.in';              % Specify the file name
ReadFile(fname);

% Write basic data of program 
WriteParasOut();

% Form the stiffness matrix
GetStiff();

% Triangularize stiffness matrix
Solve();

% Finalize
Finalize();

% ----------------------- Functions -----------------------------------

% Functions
% Add paths of functions
function AddPath()
    clear;
    close all;
    clc;
    
    addpath .\SRC\Initiation
    addpath .\SRC\BasicData
    addpath .\SRC\Mechanics
    addpath .\SRC\Mechanics\Truss
    addpath .\SRC\Mechanics\Q4
    addpath .\SRC\Mechanics\Beam
    
    addpath .\SRC\Solver
end

function Finalize()
    global cdata;
    TIM = cdata.TIM;
    time = zeros(5, 1, 'double');
    time(1) = etime(TIM(2,:), TIM(1,:));
    time(2) = etime(TIM(3,:), TIM(2,:));
    time(3) = etime(TIM(4,:), TIM(3,:));
    time(4) = etime(TIM(5,:), TIM(4,:));
    time(5) = etime(TIM(5,:), TIM(1,:));
    
    fprintf(cdata.IOUT, ['\n\n' ...
        ' S O L U T I O N   T I M E   L O G   I N   S E C\n\n' ...
        '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n' ...
        '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n' ...
        '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . = %12.2f\n' ...
        '     TIME FOR LOAD CASE SOLUTIONS  . . . . . . . . . . = %12.2f\n\n' ...
        '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n'], ...
        time(1), time(2), time(3), time(4),time(5));
    
    fprintf(['\n' ...
        ' S O L U T I O N   T I M E   L O G   I N   S E C\n\n' ...
        '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n' ...
        '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n' ...
        '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . = %12.2f\n' ...
        '     TIME FOR LOAD CASE SOLUTIONS  . . . . . . . . . . = %12.2f\n\n' ...
        '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n'], ...
        time(1), time(2), time(3), time(4),time(5));
    
    fclose(cdata.IIN);
    fclose(cdata.IOUT);
    fclose(cdata.TOUT);%(改动)
    fclose(cdata.POUT);%(改动)
end