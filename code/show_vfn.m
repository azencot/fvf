function show_vfn(mesh,vf)
if size(vf,2) == 1
    vf = reshape(vf,mesh.nf,3);
end
show_vf(mesh,vf,MESH.normv(vf));