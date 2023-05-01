%* *****************************************************************
%* - Basic data class of STAPMAT                                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Storing variables which control the running of STAPMAT      *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%* *****************************************************************
classdef ControlData
    properties
        NUMNP;         % Total number of nodal points
                       % = 0 : Program stop

        NPAR;          % Element group control data
                       %   NPAR(1) - Element type
                       %             1 : Truss element
                       %             2 /3: Q4 element                       
                       %   NPAR(2) - Number of elements
                       %   NPAR(3) - Number of different sets of material
                       %             and cross-sectional constants
        NUMEG;         % Total number of element groups, > 0
        NLCASE;        % Number of load case (>0)%载荷工况数（分布力、集中力）
        LL;            % Load case number
        NLOAD;         % The number of concentrated loads applied in this load case
        
        %Time integral ControlData
        MTIME;         % The time span of motion ; unit : s ; 0 - Static solution , Other - Time integral
        MDELTAT;       % The time step of integration; unit : s; MDELTAT(NLCASE,1)
        ALPHA;         % Rayleigh Damp parameter alpha, ALPHA(NLCASE,1)
        BETA;          % Rayleigh Damp parameter beta, BETA(NLCASE,1)
        NVEL;          % Number of initial velocity vector, NVEO(NLCASE,1)

        MODEX;         % Solution mode: 0 - data check only; 1 - execution;

        TIM;           % Timing information %求解器运行的时间
        HED;           % Master heading information for user in labeling the output %对问题的描述，无所谓，但是要有

        IIN;           % file pointer used for input 输入文件指针
        IOUT;          % file pointer used for output 输出文件指针
    end
end
