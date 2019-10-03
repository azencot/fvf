function vf_fquiver(mesh,vf,k)

if nargin < 3
    k = mesh.nf;
end

v = mesh.vertices; t = mesh.triangles;
vm = (v(t(:,1),:)+v(t(:,2),:)+v(t(:,3),:))/3;
vi = 1:mesh.nf/k:mesh.nf;


quiver3(vm(vi,1),vm(vi,2),vm(vi,3),vf(vi,1),vf(vi,2),vf(vi,3),'Color','k');
axis equal; axis off; cameratoolbar; cameratoolbar('SetCoordSys','none');