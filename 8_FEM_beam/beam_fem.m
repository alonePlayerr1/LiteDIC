% Transient FE analysis of a cantilever beam (PDE Toolbox).
% REQUIRES the Partial Differential Equation Toolbox.
%
% SPEED: the structural transient solver runs on the CPU; run time is set
% by the number of mesh elements and time steps. Coarse mesh (Hmax) and few
% time steps to keep it fast. For higher accuracy decrease Hmax (e.g. 0.002) and increase the
% number of steps in tsim. The target quantity (first natural frequency,
% ~420 Hz, matching the analytical cantilever value) converges already on
% the coarse mesh.
gm = multicuboid(0.1,0.005,0.005);
E = 210E9;
nu = 0.3;
rho = 7800;
TipVertex = addVertex(gm,Coordinates=[0.05,0,0.005]);
sModel = createpde(structural="transient-solid");
sModel.Geometry = gm;
msh = generateMesh(sModel, "Hmax", 0.004);   % coarse mesh
figure
pdegplot(sModel,FaceLabels="on");
pdemesh(msh);
title("Beam model")
drawnow
structuralProperties(sModel,YoungsModulus=E,PoissonsRatio=nu,MassDensity=rho);
structuralBC(sModel,Face=5,Constraint="fixed");
firstNF = 2639;
Tfundamental = 2*pi/firstNF;
Te = 0.02*Tfundamental;
structuralBoundaryLoad(sModel,Vertex=TipVertex,Force=[0;0;-100], ...
                              EndTime=Te,Label="force");
structuralIC(sModel,Velocity=[0;0;0]);
ncycles = 8;                                          % fewer cycles = faster
tsim = linspace(0,ncycles*Tfundamental,20*ncycles);   % fewer steps = faster
R1 = solve(sModel,tsim);
figure
plot(tsim,R1.Displacement.uz(TipVertex,:))
title('Vertical Displacement of Beam Tip')
legend('Structural PDE model')
xlabel('Time')
ylabel('Displacement')
