function PostProcessor(Node, NGauss, StressCollection)
% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
NUMNP = cdata.NUMNP;
NUME = sdata.NUME;

fprintf(IOUT, ['\n\n  S T R E S S  R E C O V E R Y  A T  N O D A L P O I N T S\n' ...
    '       NODE     Sigma_XX       Sigma_YY       Sigma_ZZ       Sigma_XY       Sigma_YZ       Sigma_ZX\n' ... 
    '       NUMBER\n']);%现在还是内力

NodeRelationFlag = zeros(NUMNP,16);
ref1 = 15;
ref2 = 16;
NStress = 6;
Stress = zeros(6,NUMNP);
for N = 1:NUME
    for i = 1:2
        j = NodeRelationFlag(Node(N,i),ref1) + 1;
        NodeRelationFlag(Node(N,i),ref1) = j;
        if j==1
            NodeRelationFlag(Node(N,i), ref2) = i;
        end
        NodeRelationFlag(Node(N,i),j) = N;
    end
end
i = max(NodeRelationFlag(:,ref1))*NGauss;
value = zeros(i,6+NStress);
for L =1:NUMNP
    Nval = NodeRelationFlag(L,ref1) * NGauss;
    if Nval~=0
        for ind2 = 1:NodeRelationFlag(L,ref1)
            N = NodeRelationFlag (L, ind2);
            ind1 = (N-1)*NGauss+1;
            for j = 1:NGauss
                ind0 = (ind2-1)*NGauss+j;
                Stress(1:NStress,L) = StressCollection (1:NStress,ind1+mod(ind0-1,NGauss));
                value(ind0,1:NStress) = Stress(1:NStress,L);
            end
        end
    end
    Stress(1:NStress,L) = (sum(value(1:Nval, 1:NStress),1))/Nval;
    
    fprintf(IOUT, ' %10d           %13.6e    %13.6e    %13.6e    %13.6e    %13.6e    %13.6e\n', L, Stress(1,L), Stress(2,L), Stress(3,L), Stress(4,L), Stress(5,L), Stress(6,L));
end
TecplotOut(Stress)
ParaviewOut(Stress)
end

function TecplotOut(Stress)

global cdata;
global sdata;

cdata.TOUT = fopen('.\Data\TECPLOT.dat', 'w');
TOUT = cdata.TOUT;
NUMNP = cdata.NUMNP;
NUME = sdata.NUME;NODEOFELE = sdata.NODEOFELE;X = sdata.X;Y = sdata.Y;Z = sdata.Z;
DIS = sdata.DIS(:,1); ID = sdata.ID;

fprintf(TOUT, 'TITLE = "Beam Element Example"\n');
fprintf(TOUT, 'VARIABLES = "X", "Y", "Z", "URX", "URY", "URZ", "FX", "FY", "FZ", "MX", "MY", "MZ"\n');
fprintf(TOUT, 'ZONE T="Beam Element", N=%5d, E=%5d, F=FEPOINT, ET=LINESEG, C=CYAN\n',NUMNP,NUME);
D = zeros(6, 1, 'double');
for i = 1:NUMNP
    D(:) = 0;
    if (ID(1, i) ~= 0) D(1) = DIS(ID(1, i)); end
    if (ID(2, i) ~= 0) D(2) = DIS(ID(2, i)); end
    if (ID(3, i) ~= 0) D(3) = DIS(ID(3, i)); end
    if (ID(4, i) ~= 0) D(4) = DIS(ID(4, i)); end
    if (ID(5, i) ~= 0) D(5) = DIS(ID(5, i)); end
    if (ID(6, i) ~= 0) D(6) = DIS(ID(6, i)); end
    fprintf(TOUT, '%18.6e%18.6e%18.6e%18.6e%18.6e%18.6e%18.6e%18.6e%18.6e%18.6e%18.6e%18.6e\n', X(i)+D(1), Y(i)+D(2), Z(i)+D(3),D(4),D(5),D(6),Stress(1,i), Stress(2,i),Stress(3,i),Stress(4,i),Stress(5,i),Stress(6,i));
end
for i = 1:NUME
    fprintf(TOUT, '%5d%5d\n', NODEOFELE(1,i), NODEOFELE(2,i));
end
end
function ParaviewOut(Stress)

global cdata;
global sdata;

cdata.POUT = fopen('.\Data\Beam3.vtk', 'w');
POUT = cdata.POUT;
NUMNP = cdata.NUMNP;
NUME = sdata.NUME;NODEOFELE = sdata.NODEOFELE;X = sdata.X;Y = sdata.Y;Z = sdata.Z;
DIS = sdata.DIS(:,1); ID = sdata.ID;

fprintf(POUT, '# vtk DataFile Version 3.0\n');
fprintf(POUT, ' BEAM \n\n');
fprintf(POUT, 'ASCII\n');
fprintf(POUT, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(POUT, 'POINTS     %5d double\n',NUMNP);
for i = 1:NUMNP
    fprintf(POUT, '%18.6e%18.6e%18.6e\n', X(i), Y(i), Z(i));
end
fprintf(POUT, 'CELLS   %5d%5d\n',NUME,3*NUME);
for i = 1:NUME
    fprintf(POUT, '2   %5d%5d\n', NODEOFELE(1,i)-1, NODEOFELE(2,i)-1);
end
fprintf(POUT, ' CELL_TYPES   %5d\n',NUME);
for i = 1:NUME
    fprintf(POUT, ' 3\n');
end
fprintf(POUT, 'CELL_DATA   %5d\n',NUME);
fprintf(POUT, 'POINT_DATA   %5d\n',NUMNP);
fprintf(POUT, 'FIELD Result          3\n');
fprintf(POUT, 'Displacement_Load_Case01            3         %5d double\n',NUMNP);
D = zeros(6, 1, 'double');
for i = 1:NUMNP
    D(:) = 0;
    if (ID(1, i) ~= 0) D(1) = DIS(ID(1, i)); end
    if (ID(2, i) ~= 0) D(2) = DIS(ID(2, i)); end
    if (ID(3, i) ~= 0) D(3) = DIS(ID(3, i)); end
    fprintf(POUT, '%18.6e%18.6e%18.6e\n', D(1), D(2), D(3));
end
fprintf(POUT, 'Rotation_Load_Case01            3       %5d double\n',NUMNP);
for i = 1:NUMNP
    D(:) = 0;
    if (ID(4, i) ~= 0) D(4) = DIS(ID(4, i)); end
    if (ID(5, i) ~= 0) D(5) = DIS(ID(5, i)); end
    if (ID(6, i) ~= 0) D(6) = DIS(ID(6, i)); end
    fprintf(POUT, '%18.6e%18.6e%18.6e\n', D(4),D(5),D(6));
end
fprintf(POUT, 'Stress_Load_Case01            6       %5d double\n',NUMNP);
for i = 1:NUMNP
    fprintf(POUT, '%18.6e%18.6e%18.6e%18.6e%18.6e%18.6e\n', Stress(1,i), Stress(2,i),Stress(3,i),Stress(4,i),Stress(5,i),Stress(6,i));
end
end