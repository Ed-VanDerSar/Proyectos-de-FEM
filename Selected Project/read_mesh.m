

function [mesh] = read_mesh(mesh_filename)
% get_faces:
% Input: Nombre del archivo .m creado usando gmsh
% Output:  Estructura mesh que tiene:
%
%          nb_nodes            : numero de nodos de la malla
%          nb_elems            : numero de triangulos/elementos de la malla
%          nodes               : nodos de la mallas en 2D
%          nb_faces            : numero de caras
%          nb_ofaces           : numero de caras en la frontera
%          faces_nodes_conn    : conectividad caras->nodos
%          ofaces_nodes_conn   : conectividad caras en la frontera->nodos
%          ofaces_bool         : vector de booleanos que indica si una cara
%                                esta en la frontera
%          onodes_bool         : vector de booleanos que indica si un node
%                                esta en la frontera
%          elems_nodes_conn    : conectividad elementos->nodos
%          elems_faces_conn    : conectividad elementos->caras
%          faces_elems_conn    : conectividad caras->elementos




run(mesh_filename); %importa  estructura msh

mesh.nb_nodes=msh.nbNod;
mesh.nodes=msh.POS;
mesh.elems_nodes_conn=msh.TRIANGLES;

foo=size(msh.TRIANGLES);
mesh.nb_elems=foo(1);
foo=size(msh.LINES);
mesh.nb_ofaces=foo(1); %cara externa

sorted_triangles=zeros(mesh.nb_elems,3);
mesh.elems_faces_conn=zeros(mesh.nb_elems,3);
mesh.faces_elems_conn=zeros(0,2);
mesh.faces_nodes_conn=zeros(0,2);
mesh.ofaces_bool=zeros(0,2);
for elem_idx=1:mesh.nb_elems 

    sorted_triangles(elem_idx,:)=sort(msh.TRIANGLES(elem_idx,1:3));
    lfaces=zeros(3,2);
    lfaces(1,:)=sorted_triangles(elem_idx,1:2);
    lfaces(2,:)=sorted_triangles(elem_idx,2:3);
    lfaces(3,:)=sorted_triangles(elem_idx,1:2:3);

    for lface_idx=1:3

    cface=lfaces(lface_idx,:);
    gfidx_v=(mesh.faces_nodes_conn(:,1)==cface(1)) & (mesh.faces_nodes_conn(:,2)==cface(2));
    gfidx=find(gfidx_v,1);
        if isempty(gfidx) %cface was not found in global faces
            mesh.faces_nodes_conn(end+1,:)=cface;
            gfaces_size=size(mesh.faces_nodes_conn);
            mesh.elems_faces_conn(elem_idx,lface_idx)=gfaces_size(1);

            mesh.faces_elems_conn(end+1,1)=elem_idx; %add new triangle to face connectivity 
        else
            mesh.elems_faces_conn(elem_idx,lface_idx)=gfidx;
            mesh.faces_elems_conn(gfidx,2)=elem_idx; %add new triangle to face connectivity
        end

    end %j

end %nb_elems
mesh.ofaces_bool=mesh.faces_elems_conn(:,2)==0;
mesh.ofaces_nodes_conn=mesh.faces_nodes_conn(mesh.ofaces_bool,:);

foo=size(mesh.faces_elems_conn);
mesh.nb_faces=foo(1);


%Build vector booleans for nodes at the boundary
mesh.onodes_bool=false(mesh.nb_nodes,1);


for i=1:mesh.nb_ofaces
    mesh.onodes_bool(mesh.ofaces_nodes_conn(i,1))=true;
    mesh.onodes_bool(mesh.ofaces_nodes_conn(i,2))=true;
end

clear msh %free tmp structure

end
