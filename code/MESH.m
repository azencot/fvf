classdef MESH < handle
    
    properties (Access='public')
        
        name
        
        vertices
        triangles
        
        nv              % #vertices
        nf              % #faces
        
        va              % vertex areas
        ta              % face areas
        
        Nf              % normal-per-face
        
        G               % tangential gradient operator
        D               % tangential divergence operator
    end
    
    properties (Access='protected')
        
        E1, E2, E3      % triangle edges

    end
    
    methods
        
        function [ mesh ] = MESH( meshname )
            
            if nargin < 1; meshname = 'sphere_s3'; end
                        
            mesh.name = meshname;
            
            [mesh.vertices, mesh.triangles] = MESH_READER.readOff([meshname '.off']);
            
            mesh.nv = size(mesh.vertices,1);
            mesh.nf = size(mesh.triangles,1);
            
            mesh.Nf = face_normals( mesh );
            mesh.ta = face_areas( mesh );
            mesh.va = vertex_areas( mesh );
            
            [mesh.E1,mesh.E2,mesh.E3] = face_edges( mesh );
                        
            mesh.G = grad( mesh );
            mesh.D = div( mesh );
        end
        
        function [ ta ] = face_areas( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            ta = mesh.normv( cross( P1, P2 ) ) / 2;
        end
        
        function [ va ] = vertex_areas( mesh )
            va = full( sum( mass_matrix(mesh), 2 ));
        end
        
        function [ M ] = mass_matrix( mesh )
            T = mesh.triangles; 
            
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            Mij = 1/12*[mesh.ta; mesh.ta; mesh.ta];
            Mji = Mij;
            Mii = 1/6*[mesh.ta; mesh.ta; mesh.ta];
            In = [I;J;I];
            Jn = [J;I;I];
            Mn = [Mij;Mji;Mii];
            M = sparse(In,Jn,Mn,mesh.nv,mesh.nv);
        end
                
        function [ Nf ] = face_normals( mesh )
            
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            Nf = cross( P1, P2 );
            Nf = MESH.normalize_vf( Nf );
        end
               
        function [ E1, E2, E3 ] = face_edges( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            E1 = X(T(:,3),:) - X(T(:,2),:);
            E2 = X(T(:,1),:) - X(T(:,3),:);
            E3 = X(T(:,2),:) - X(T(:,1),:);
        end
        
        function [ G ] = grad( mesh )
            % G corresponds to eq. (3.9) in Polygon mesh processing book
            I = repmat(1:mesh.nf,3,1);
            II = [I(:); I(:)+mesh.nf; I(:)+2*mesh.nf];

            J = mesh.triangles';
            JJ = [J(:); J(:); J(:)];

            RE1 = mesh.rotate_vf( mesh.E1 );
            RE2 = mesh.rotate_vf( mesh.E2 );
            RE3 = mesh.rotate_vf( mesh.E3 );

            S = [ RE1(:) RE2(:) RE3(:) ]';
            SS = S(:);

            G = sparse(II,JJ,SS,3*mesh.nf,mesh.nv);

            TA = .5 * repmat(1 ./ mesh.ta,1,3);
            A = spdiags(TA(:),0,3*mesh.nf,3*mesh.nf);

            G = A*G;
        end
        
        function [ D ] = div( mesh )
            % D corresponds to eq. (3.12) in Polygon mesh processing book
            IVA = spdiags(1./mesh.va,0,mesh.nv,mesh.nv);
            TAC = spdiags(repmat(mesh.ta,3,1),0,3*mesh.nf,3*mesh.nf);

            D = - IVA * mesh.G' * TAC;
        end
        
        function [ op ] = vf2op( mesh, vf )
            vf = reshape( vf, mesh.nf, 3 );
            
            RE1 = mesh.rotate_vf( mesh.E1 );
            RE2 = mesh.rotate_vf( mesh.E2 );
            RE3 = mesh.rotate_vf( mesh.E3 );
            
            T = mesh.triangles;
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            Sij = 1/6*[dot(RE2,vf,2); dot(RE3,vf,2); dot(RE1,vf,2)];
            Sji = 1/6*[dot(RE1,vf,2); dot(RE2,vf,2); dot(RE3,vf,2)];
            In = [I;J;I;J];
            Jn = [J;I;I;J];
            Sn = [Sij;Sji;-Sij;-Sji];
            W = sparse(In,Jn,Sn,mesh.nv,mesh.nv);
            IVA = spdiags(1./mesh.va,0,mesh.nv,mesh.nv);
            
            op = IVA*W;
        end
               
        function [ vf ] = project_vf( mesh, vf )
            vf = vf - repmat(dot(vf,mesh.Nf,2),1,3).*mesh.Nf;
        end
        
        function [ rvf ] = rotate_vf( mesh, vf )
            rvf = cross( mesh.Nf, vf );
        end
        
    end
    
    methods (Static)
        
        function [ nv ] = normv( v )
            nv = sqrt(sum(v.^2,2));
        end
        
        function [ nnv ] = normalize_vf( v )
            nnv = v ./ repmat(MESH.normv(v),1,3);
        end

    end
end