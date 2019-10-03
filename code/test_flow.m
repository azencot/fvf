%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

paths

t = 2; k = 100; h = t/k;

meshname = 'sphere_s2';
mesh = MESH(meshname);

run_comp = 1;
run_vis = 1;

%%%%%%%%%%%%%%%
% Computation %
%%%%%%%%%%%%%%%
if run_comp == 1

% eigenfunctions of the Laplacian are smooth and nice for transport
[evecs,evals] = eigs(mesh.D*mesh.G,9,1e-5);
f = evecs(:,7);

vf = mesh.project_vf( repmat( [0 1 0], mesh.nf, 1 ) );
op = mesh.vf2op( vf );

F = zeros(mesh.nv,k); F(:,1) = f;
for i = 2:k
    F(:,i) = expmv(-h,op,F(:,i-1),[],'double');
end

end

%%%%%%%%%%%%%%%%%
% Visualization %
%%%%%%%%%%%%%%%%%
if run_vis == 1
    
figure; show_vfn(mesh,vf);
figure;    
for i = 1:k
    clf; show_func(mesh,F(:,i));
    pause(0.01);    
end

end