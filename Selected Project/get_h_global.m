% get_h_global: regresa la h global de la malla
% Input: estructura msh 
% Output: h_global


function [h_global] = get_h_global(msh)

% guardamos/renombramos variables importantes

h_global=0.0;

for i=1:msh.nb_elems 

    %diferencia entre los puntos de cada cara
    dface_1=msh.nodes(msh.elems_nodes_conn(i,1),1:2)-msh.nodes(msh.elems_nodes_conn(i,2),1:2);
    dface_2=msh.nodes(msh.elems_nodes_conn(i,2),1:2)-msh.nodes(msh.elems_nodes_conn(i,3),1:2);
    dface_3=msh.nodes(msh.elems_nodes_conn(i,3),1:2)-msh.nodes(msh.elems_nodes_conn(i,1),1:2);
    h_T=max(max(norm(dface_1),norm(dface_2)),norm(dface_3));
    h_global=max(h_global,h_T);
end

end
