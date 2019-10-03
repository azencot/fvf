function show_vf(mesh,vf,f,bn)

hold on;
if nargin < 3
    patch('faces',mesh.triangles,'vertices',mesh.vertices,...
        'FaceColor',[1,0.936,0.859],'EdgeColor',[0.5,0.5,0.5]);
else
    show_func(mesh,f);
end

if nargin < 4
    bn = 0;
end

if bn == 1
    vf = normalize_vf(vf);
end

vf_fquiver(mesh,vf,mesh.nf);
colorbar; hold off;