%* *****************************************************************
%* - Basic data class of STAPMAT                                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Storing variables used in solving process                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%* *****************************************************************

classdef SolutionData
    properties (Constant)
        % Gauss coord, 1D to 3D
        GC1 = double(0.0);
        GC2 = double([1/3,-1/3]);
        GC3 = double([sqrt(0.6), 0.0, -sqrt(0.6)]);
        % Gauss weight, 1D to 3D
        GW1 = double(2.0);
        GW2 = double([1.0, 1.0]);
        GW3 = double([5.0/9.0, 8.0/9.0, 5.0/9.0]);
    end
    properties
        % Basic data
        ID;       % int, ID(3, NUMNP), Boundary condition codes (0=free, 1=deleted)
        IDOrigin; % int, backups of ID after computing of NEQ
        X;        % double, X(NUMNP), X coordinates
        Y;        % double, Y(NUMNP), Y coordinates
        Z;        % double, Z(NUMNP), Z coordinates
        R;        % double, R(NEQ), Load vector
        Q;        % double, Q(NEQ), Back up of Load vector
        NOD;      % int, NOD(NLOAD), Node number to which this load is applied (1~NUMNP)
        IDIRN;    % int, IDIRN(NLOAD), Degree of freedom number for this load component
                  %                     1 : X-direction;
                  %                     2 : Y-direction;
                  %                     3 : Z-direction;
        FLOAD;    % double, FLOAD(NLOAD), Magnitude of load
           
        
        
        % Element data
        NUME;     % int, number of elements    %单元个数
        NNODE;    % int, number of nodes in an element
        NINIP;    % int, number of integration points in an element
        NDOF;     % int, the DOF of displacement
        NODEOFELE;%（改动）
        NSTIFF;   % int, the number of number in element stiffness matrix
        XYZ;      % double, XYZ(NDOF*NNODE, NUME), element position
        POS;      % double ,POS(NEQ) ,free node position at T= 0;%(改动)
        
        InitCoord;  % double array, integration coordinates
        InitWeight; % double array, integration weights
        
        % Material data
        NUMMAT;     % int, the number of types of material  %不同的材料种类
        E;          % double array, Young's Modulus
        AREA;       % double array, Cross sectional area of truss
        Iy;         %截面惯性矩Iy%(改动)
        Iz;         %截面惯性矩Iz
        Jx;         %截面极惯性矩Jx
        NU;         % double array, Possion ratio
        RHO;        % double array, Density
        MU;         % double array, Damp
        MATP;       % int, MATP(NUME), types of elements
        
        % Solve data
        NEQ;        % int, Number of equations
        NWK;        % Number of matrix elements
        MK;         % Maximum half bandwidth
        MHT;        % int, MHT(NEQ), Vector of column heights
        LM;         % int, LM(NDOF*NNODE, NUME), Connectivity matrix
        MAXA;       % int, MAXA(NEQ)
        STIFF;      % double, STIFF(NWK), store the elements of stiffness matrix
        STIFFOrigin;% double, sparse,STIFF(NWK), Back up of the Origin Stiff%(改动)
        MASS = 0;   % double, sparse,MASS(NEQ), store the elements of initial mass matrix%(改动)
        DAMP = 0;   % double, sparse,DAMP(NEQ), store the elements of initial damp matrix%(改动)
        
        %Initial conditions%(改动)
        NODVEL;   % int, NODVEL(NVEL), Node number to which this velocity is applied (1~NUMNP)
        IDVEL;    % int, IDVEL(NVEL), Degree of freedom number for this velocity component
                  %                     1 : X-direction;
                  %                     2 : Y-direction;
                  %                     3 : Z-direction;
        MVEL;     % double, DVEL(NVEL), Magnitude of velocity
        V;        % double, V(NEQ), velocity vector

        % Result data
        DIS;      % double, DIS(NEQ, NLCASE), Displacement of nodes
        STRAIN;   % double, STRAIN(NEQ, NLCASE), Strain
        STRESS;   % double, STRESS(NEQ, NLCASE), Stress
        
    end
end