%* *****************************************************************
%* - Function of STAPMAT in initialization phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read input file of STAPMAT                                  *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Initiation/ReadFile.m - InitBasicData()                 *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function ReadFile(fname)
fname = strcat('.\Data\', fname);           % Deal the filename

% Get global class
global cdata;
global sdata;

% Open files
cdata.IIN = fopen(fname, 'r');

% Begin Read input file
fprintf('Input phase ...\n\n');

% the first time stamp
cdata.TIM = zeros(5, 6, 'double');
cdata.TIM(1,:) = clock;

IIN = cdata.IIN;
%% Read Control data
cdata.HED = fgetl(IIN);

tmp = str2num(fgetl(IIN));
cdata.NUMNP = tmp(1);
cdata.NUMEG = tmp(2);
cdata.NLCASE = tmp(3);
cdata.MODEX = tmp(4);

if (cdata.NUMNP == 0) return; end

%% Read nodal point data(这里有修改)
InitBasicData();
% Define local variables to speed
ID = sdata.ID; X = sdata.X; Y = sdata.Y; Z = sdata.Z;
for i = 1:cdata.NUMNP
    tmp = str2num(fgetl(IIN));
    ID(1, i) = tmp(2);
    ID(2, i) = tmp(3);
    ID(3, i) = tmp(4);
    ID(4, i) = tmp(5);
    ID(5, i) = tmp(6);
    ID(6, i) = tmp(7);
    X(i) = tmp(8);
    Y(i) = tmp(9);
    Z(i) = tmp(10);
end
sdata.ID = ID; sdata.X = X; sdata.Y = Y; sdata.Z = Z;
%% Compute the number of equations(这里有修改)
sdata.IDOrigin = ID;
NEQ = 0;
for N=1:cdata.NUMNP
    for I=1:6
        if (ID(I,N) == 0)
            NEQ = NEQ + 1;
            ID(I,N) = NEQ;
        else
            ID(I,N) = 0;
        end
    end
end
sdata.ID = ID;
sdata.NEQ = NEQ;
%% Read load data(这里有修改)
% Init control data
NLCASE = cdata.NLCASE;
sdata.R = zeros(NEQ, NLCASE);
R = sdata.R;
% Read data
for N = 1:cdata.NLCASE
    tmp = str2num(fgetl(IIN));
    cdata.LL = tmp(1); cdata.NLOAD = tmp(2);
    NLOAD = cdata.NLOAD;
%   Init load data
    sdata.NOD = zeros(NLOAD, 1);
    sdata.IDIRN = zeros(NLOAD, 1);
    sdata.FLOAD = zeros(NLOAD, 1);
    NOD = sdata.NOD; IDIRN = sdata.IDIRN; FLOAD = sdata.FLOAD;
    
%   Read load data
    for I = 1:NLOAD
        tmp = str2num(fgetl(IIN));
        NOD(I) = tmp(1);
        IDIRN(I) = tmp(2);
        FLOAD(I) = tmp(3);
    end
    if (cdata.MODEX == 0) return; end
    
%   Compute load vector
    for L = 1:NLOAD
        II = ID(IDIRN(L), NOD(L));
        if (II > 0) R(II, N) = R(II, N) + FLOAD(L); end
    end
    sdata.NOD = NOD; sdata.IDIRN = IDIRN; sdata.FLOAD = FLOAD; sdata.R = R;
end

end

%% Functions


% InitBasicData
function InitBasicData()
global cdata;
global sdata;

cdata.NPAR = zeros(10, 1);

sdata.ID = zeros(3,cdata.NUMNP);
sdata.X = zeros(cdata.NUMNP, 1);
sdata.Y = zeros(cdata.NUMNP, 1);
sdata.Z = zeros(cdata.NUMNP, 1);
end