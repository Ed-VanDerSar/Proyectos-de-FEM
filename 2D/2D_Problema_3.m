%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Teoría y Práctica Elementos Finitos
% Posgrado en Matematicas-UNAM-CdMx
% Presenta: Luis Edwin Aguilar Anzures
% email: eda.math10@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descripcion: Dada u calculamos el error en L^2 y H^1 y sus respectivas
% tasas de decrecimiento para las mallas dadas entre ella
% misma y la solución aproximada obtenida con el elemento finito de
% Crouseix-Raviart de orden 1 en (0,1)x(0,1).
% Ouput: Tabla que contiene el error en norma L^2, norma H^1 y sus
% respectivas tasas de convergencia para una función dada y su interpolador
% respecto a los elementos finitos de Crouseix-Raviart de primer orden en (0,1)x(0,1) en las mallas
% dadas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cargamos las mallas dadas
msh_filenames = ["square_1.m" "square_2.m" "square_3.m"...
  "square_4.m" "square_5.m" "square_6.m"];
h_vec=zeros(1,6); % Vector que contiene el h para cada malla
L2_error_vec=zeros(1,6); % Vector que contiene el error en L2 para cada malla
H1_error_vec=zeros(1,6); % Vector que contiene el error en H1 para cada ciclo
u_x_y=@u_exact; %cargamos la función dada  u


%empieza el ciclo para cada malla dada.
for j=1:6
msh=read_mesh(msh_filenames(j));
h_vec(j)=get_h_global(msh);
dofs=u_x_y(msh.nodes(:,1),msh.nodes(:,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculo del eror u - u_h en la norma L2

%sumamos los errores L^2 en cada elelemnto para obtener el error total de
%la malla 

    L2_error_f= @L2_error; %función que calcula el error en cada elemento. Vea descripción de la función.
    error=0;
    %Ciclo para cad amalla
    for i=1:msh.nb_elems  
       error=error+ L2_error_f(i,msh);
    end
    error= sqrt(error);
    L2_error_vec(j)=error; %Guardar el valor del error L2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Calculo del eror u - u_h en la norma H1

%sumamos los errores de la seminorma de siobolev H^1 en cada elemento para obtener el error total de
%la malla y al final le sumamos el respectivo error en L^2 para calcular el
%error en la norma de Sobolev H^1.

H1_error_f= @semH1_error ; %función que calcula el error en seminorma en cada elemento. Vea descripción de la función.
   h1error=0;
    %Ciclo para cada subintervalo
    for i=1:msh.nb_elems  
        h1error= h1error+ H1_error_f(i,msh);  %suma el valor de la integral en cada subintervalo 
    end
    h1error= sqrt(h1error+(L2_error_vec(j))^2);
    H1_error_vec(j)=h1error; %Guardar el valor del error H1
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Crear tabla con valores: N| error L2| tasa_de_decrecimiento L2| error H1|
% tasa_de_decrecimiento H1

err_rate=zeros(1,6); %tasa de decrecimiento del error
err_rate(1)=1;
for i=2:6
    err_rate(i)=log(L2_error_vec(i)./L2_error_vec(i-1))./log(1/2);
end

% decrecimiento del error H1
H1err_rate=zeros(1,6); %tasa de decrecimiento del error
H1err_rate(1)=1;
for i=2:6
    H1err_rate(i)=log(H1_error_vec(i)./H1_error_vec(i-1))./log(1/2);
end

%renombramos lo títulos de la tabla: 
h_vec=h_vec';
L2_err_norm=L2_error_vec';
L2_err_rate= err_rate';
H1_err_norm=H1_error_vec';
H1_err_rate= H1err_rate';


output_table=table(h_vec,L2_err_norm,L2_err_rate,H1_err_norm,H1_err_rate) ; % Crea tabla con vectores columna                                           
disp(output_table);%Imprime en la terminal



%%%%%%%%%%%%%% FUNCIONES %%%%%%%%%%%%%%%%%

%la sol teórica 
function [val] = u_exact(x,y)
    val = cos(x*pi*4).*(cos(y*pi*4)).^2; 
end
%la derivada parcial respecto a x sol teórica
function [val] = dx_u_exact(x,y)
    val = -4*pi*sin(4*pi*x)*(cos(y*pi*4))^2; 
end
%la derivada parcial respecto a y sol teórica
function [val] = dy_u_exact(x,y)
    val = -8*pi*cos(4*pi*x)*cos(y*pi*4)*sin(y*pi*4); 
end



%%DESCRIPCIÓN: Calcularmos el interpolador en un elemento i usando las funciones de
%%forma del elemento de referencia, una transformación al elemento de
%%referencia y cuadratura gaussiana para calcular los grados de libertad
%%que dependen de integrales en las caras.
function [val] = eval_hat_functions(x,y,i,msh)
V_1=msh.nodes(msh.elems_nodes_conn(i,1),1:2); %primer vértice del elemento
V_2=msh.nodes(msh.elems_nodes_conn(i,2),1:2); %segundo vértice del elemento
V_3=msh.nodes(msh.elems_nodes_conn(i,3),1:2); %tercer vértice del elemento
IN=inpolygon(x,y,[V_1(1),V_2(1),V_3(1)],[V_1(2),V_2(2),V_3(2)]) ; %revisa que el punto (x,y) yace en el elemento
A=[V_3'-V_1',V_2'-V_1']; %matriz de la tranfimación afin del elemento de referencia al elemento i
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)]; %matriz de la tranfimación afin del elemento i al elemento de referencia 
%puntos sobre las caras para hacer la integracuión por cuadratura gaussiana
%en cada una 
F_1_1=V_2+(V_3-V_2)*(0.5-0.5*sqrt(3/5)); %primer punto en F_1
F_1_2=V_2+(V_3-V_2)*(0.5); %segundo  punto en F_1
F_1_3=V_2+(V_3-V_2)*(0.5+0.5*sqrt(3/5)); %tercer punto en F_1
F_2_1=V_3+(V_1-V_3)*(0.5-0.5*sqrt(3/5));
F_2_2=V_3+(V_1-V_3)*(0.5);
F_2_3=V_3+(V_1-V_3)*(0.5+0.5*sqrt(3/5));
F_3_1=V_1+(V_2-V_1)*(0.5-0.5*sqrt(3/5));
F_3_2=V_1+(V_2-V_1)*(0.5);
F_3_3=V_1+(V_2-V_1)*(0.5+0.5*sqrt(3/5));
%integrales en las caras usando cuadratura Gaussiana a tres puntos en el intervalo [0,1],
%(con el cual se parametriza la cara) tales integrales son los grados de
%libertad asociados a las caras de este elemento.
IntF_1= [5/18,4/9,5/18]*[u_exact(F_1_1(1),F_1_1(2)),u_exact(F_1_2(1),F_1_2(2)),u_exact(F_1_3(1),F_1_3(2))]';
IntF_2= [5/18,4/9,5/18]*[u_exact(F_2_1(1),F_2_1(2)),u_exact(F_2_2(1),F_1_2(2)),u_exact(F_2_3(1),F_2_3(2))]';
IntF_3= [5/18,4/9,5/18]*[u_exact(F_3_1(1),F_3_1(2)),u_exact(F_3_2(1),F_3_2(2)),u_exact(F_3_3(1),F_3_3(2))]';
%coordenadas baricentricas para hacer la funci{on de forma
%\theta_i=1-2*\lambda_i
 lambda1=(1-(x-V_1(1))*A_in(1,1)-(y-V_1(2))*A_in(1,2)-(x-V_1(1))*A_in(2,1)-(y-V_1(2))*A_in(2,2) );
 lambda2=((x-V_1(1))*A_in(1,1)+(y-V_1(2))*A_in(1,2));
 lambda3=((x-V_1(1))*A_in(2,1)+(y-V_1(2))*A_in(2,2));
 %interpolador I=\Sum_[i=1:3]IntF_i\theta_i
 val=IntF_1*(1-2*lambda1)+IntF_2*(1-2*lambda2)+IntF_3*(1-2*lambda3) ;
val=IN.*val;
end



%%DESCRIPCIÓN: Calcularmos la derivada parcial respecto a la coordenada x interpolador(x,y) de en un elemento i usando las funciones de
%%forma del elemento de referencia y una transformación al elemento de
%%referencia 
function [val] = dx_eval_hat_functions(x,y,i,msh)

V_1=msh.nodes(msh.elems_nodes_conn(i,1),1:2);
V_2=msh.nodes(msh.elems_nodes_conn(i,2),1:2);
V_3=msh.nodes(msh.elems_nodes_conn(i,3),1:2);
IN=inpolygon(x,y,[V_1(1),V_2(1),V_3(1)],[V_1(2),V_2(2),V_3(2)]) ;
A=[V_3'-V_1',V_2'-V_1'];
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)];
%De nuevo necesitamos calcular los grados de libertad:
%puntos sobre las caras para hacer la integracuión por cuadratura gaussiana
%en cada una 
F_1_1=V_2+(V_3-V_2)*(0.5-0.5*sqrt(3/5)); %primer punto en F_1
F_1_2=V_2+(V_3-V_2)*(0.5); %segundo  punto en F_1
F_1_3=V_2+(V_3-V_2)*(0.5+0.5*sqrt(3/5)); %tercer punto en F_1
F_2_1=V_3+(V_1-V_3)*(0.5-0.5*sqrt(3/5));
F_2_2=V_3+(V_1-V_3)*(0.5);
F_2_3=V_3+(V_1-V_3)*(0.5+0.5*sqrt(3/5));
F_3_1=V_1+(V_2-V_1)*(0.5-0.5*sqrt(3/5));
F_3_2=V_1+(V_2-V_1)*(0.5);
F_3_3=V_1+(V_2-V_1)*(0.5+0.5*sqrt(3/5));
%integrales en las caras usando cuadratura Gaussiana a tres puntos en el intervalo [0,1],
%(con el cual se parametriza la cara) tales integrales son los grados de
%libertad asociados a las caras de este elemento.
IntF_1= [5/18,4/9,5/18]*[u_exact(F_1_1(1),F_1_1(2)),u_exact(F_1_2(1),F_1_2(2)),u_exact(F_1_3(1),F_1_3(2))]';
IntF_2= [5/18,4/9,5/18]*[u_exact(F_2_1(1),F_2_1(2)),u_exact(F_2_2(1),F_1_2(2)),u_exact(F_2_3(1),F_2_3(2))]';
IntF_3= [5/18,4/9,5/18]*[u_exact(F_3_1(1),F_3_1(2)),u_exact(F_3_2(1),F_3_2(2)),u_exact(F_3_3(1),F_3_3(2))]';
%directamente calculamos la parcial \partial x  del interpolador
val=-2*IntF_1*(-A_in(1,1)-A_in(2,1));
val = val -2*IntF_2*A_in(1,1);
val = val -2*IntF_3*A_in(2,1);
val=IN.*val;
end

%%DESCRIPCIÓN: Calcularmos la derivada parcial respecto a la coordenada y interpolador(x,y) de en un elemento i usando las funciones de
%%forma del elemento de referencia y una transformación al elemento de
%%referencia 
function [val] = dy_eval_hat_functions(x,y,i,msh)
V_1=msh.nodes(msh.elems_nodes_conn(i,1),1:2);
V_2=msh.nodes(msh.elems_nodes_conn(i,2),1:2);
V_3=msh.nodes(msh.elems_nodes_conn(i,3),1:2);
IN=inpolygon(x,y,[V_1(1),V_2(1),V_3(1)],[V_1(2),V_2(2),V_3(2)]) ;
A=[V_3'-V_1',V_2'-V_1'];
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)];
%De nuevo necesitamos calcular los grados de libertad:
%puntos sobre las caras para hacer la integracuión por cuadratura gaussiana
%en cada una 
F_1_1=V_2+(V_3-V_2)*(0.5-0.5*sqrt(3/5)); %primer punto en F_1
F_1_2=V_2+(V_3-V_2)*(0.5); %segundo  punto en F_1
F_1_3=V_2+(V_3-V_2)*(0.5+0.5*sqrt(3/5)); %tercer punto en F_1
F_2_1=V_3+(V_1-V_3)*(0.5-0.5*sqrt(3/5));
F_2_2=V_3+(V_1-V_3)*(0.5);
F_2_3=V_3+(V_1-V_3)*(0.5+0.5*sqrt(3/5));
F_3_1=V_1+(V_2-V_1)*(0.5-0.5*sqrt(3/5));
F_3_2=V_1+(V_2-V_1)*(0.5);
F_3_3=V_1+(V_2-V_1)*(0.5+0.5*sqrt(3/5));
%integrales en las caras usando cuadratura Gaussiana a tres puntos en el intervalo [0,1],
%(con el cual se parametriza la cara) tales integrales son los grados de
%libertad asociados a las caras de este elemento.
IntF_1= [5/18,4/9,5/18]*[u_exact(F_1_1(1),F_1_1(2)),u_exact(F_1_2(1),F_1_2(2)),u_exact(F_1_3(1),F_1_3(2))]';
IntF_2= [5/18,4/9,5/18]*[u_exact(F_2_1(1),F_2_1(2)),u_exact(F_2_2(1),F_1_2(2)),u_exact(F_2_3(1),F_2_3(2))]';
IntF_3= [5/18,4/9,5/18]*[u_exact(F_3_1(1),F_3_1(2)),u_exact(F_3_2(1),F_3_2(2)),u_exact(F_3_3(1),F_3_3(2))]';
%directamente calculamos la parcial \partial y  del interpolador 
val=-2*IntF_1*(-A_in(1,2)-A_in(2,2));
val = val -2*IntF_2*A_in(1,2);
val = val -2*IntF_3*A_in(2,2);
val=IN.*val;
end

%DESCRIPCIÓN: función que mapea un punto en el elementos de referencia al elemento i
%para hacer la integración T(x,y)=A(x,y)+V_1
function [val] =map(x,y,i,msh)
%msh=read_mesh(['square_',num2str(j),'.m']);
val=[x*(msh.nodes(msh.elems_nodes_conn(i,3),1)-msh.nodes(msh.elems_nodes_conn(i,1),1))+y*(msh.nodes(msh.elems_nodes_conn(i,2),1)-msh.nodes(msh.elems_nodes_conn(i,1),1))+msh.nodes(msh.elems_nodes_conn(i,1),1),x*(msh.nodes(msh.elems_nodes_conn(i,3),2)-msh.nodes(msh.elems_nodes_conn(i,1),2))+y*(msh.nodes(msh.elems_nodes_conn(i,2),2)-msh.nodes(msh.elems_nodes_conn(i,1),2))+msh.nodes(msh.elems_nodes_conn(i,1),2) ];
end


 %DESCRIPCIÓN:Calcula el error en la seminorma de Sobolvev en el elemnto i usando un cambio de variable al elemento 
 % de referencia y el método de cuadratura Gaussiana de orden dos (tres puntos)
function [val]=semH1_error(i,msh)
S =msh.nodes(msh.elems_nodes_conn(i,1),1)*(msh.nodes(msh.elems_nodes_conn(i,2),2)-msh.nodes(msh.elems_nodes_conn(i,3),2))+msh.nodes(msh.elems_nodes_conn(i,2),1)*(msh.nodes(msh.elems_nodes_conn(i,3),2)-msh.nodes(msh.elems_nodes_conn(i,1),2))+msh.nodes(msh.elems_nodes_conn(i,3),1)*(msh.nodes(msh.elems_nodes_conn(i,1),2)-msh.nodes(msh.elems_nodes_conn(i,2),2));
 S=0.5*S; %área del elemnto i 
X_1=map(0,0.5,i,msh); %primer punto para la cuadratura
X_2=map(0.5,0.5,i,msh);  %segundo punto para la cuadratura
X_3=map(0.5,0,i,msh);%tercer punto para la cuadratura
%determinante de la derivada del mapeo que va del simplejo de referencia al
%elekemnto i 
detJ=det([msh.nodes(msh.elems_nodes_conn(i,3),1:2)'-msh.nodes(msh.elems_nodes_conn(i,1),1:2)',msh.nodes(msh.elems_nodes_conn(i,2),1:2)'-msh.nodes(msh.elems_nodes_conn(i,1),1:2)'] ) ;
%evaluamos la integral de acuerdo con la fórmula de la cuadratura
%gaussiana
val= (dx_u_exact(X_1(1),X_1(2))- dx_eval_hat_functions(X_1(1),X_1(2),i,msh)).^2+(dy_u_exact(X_1(1),X_1(2))- dy_eval_hat_functions(X_1(1),X_1(2),i,msh)).^2;
  val= val + (dx_u_exact(X_2(1),X_2(2))- dx_eval_hat_functions(X_2(1),X_2(2),i,msh)).^2+(dy_u_exact(X_2(1),X_2(2))- dy_eval_hat_functions(X_2(1),X_2(2),i,msh)).^2;
  val= val + (dx_u_exact(X_3(1),X_2(2))- dx_eval_hat_functions(X_3(1),X_3(2),i,msh)).^2+(dy_u_exact(X_3(1),X_2(2))- dy_eval_hat_functions(X_3(1),X_3(2),i,msh)).^2;
val=abs(detJ)*S*(1/3)*val;
end



%DESCRIPCIÓN:Calculla el error en la seminorma de Sobolvev en el elemnto i usando un cambio de variable al elemento 
 % de referencia y el método de cuadratura Gaussiana de orden dos (tres puntos)
function [val] =L2_error(i,msh)

 S =msh.nodes(msh.elems_nodes_conn(i,1),1)*(msh.nodes(msh.elems_nodes_conn(i,2),2)-msh.nodes(msh.elems_nodes_conn(i,3),2))+msh.nodes(msh.elems_nodes_conn(i,2),1)*(msh.nodes(msh.elems_nodes_conn(i,3),2)-msh.nodes(msh.elems_nodes_conn(i,1),2))+msh.nodes(msh.elems_nodes_conn(i,3),1)*(msh.nodes(msh.elems_nodes_conn(i,1),2)-msh.nodes(msh.elems_nodes_conn(i,2),2));
 S=0.5*S;  %área del elemnto i 
X_1=map(0,0.5,i,msh); %primer punto para la cuadratura
X_2=map(0.5,0.5,i,msh);  %segundo punto para la cuadratura
X_3=map(0.5,0,i,msh);%tercer punto para la cuadratura
%determinante de la derivada del mapeo que va del simplejo de referencia al
%elekemnto i 
detJ=det([msh.nodes(msh.elems_nodes_conn(i,3),1:2)'-msh.nodes(msh.elems_nodes_conn(i,1),1:2)',msh.nodes(msh.elems_nodes_conn(i,2),1:2)'-msh.nodes(msh.elems_nodes_conn(i,1),1:2)'] ) ;
%evaluamos la integral de acuerdo con la fórmula de la cuadratura
%gaussiana
val= (u_exact(X_1(1),X_1(2))- eval_hat_functions(X_1(1),X_1(2),i,msh)).^2;
  val= val + (u_exact(X_2(1),X_2(2))- eval_hat_functions(X_2(1),X_2(2),i,msh)).^2;
  val= val + (u_exact(X_3(1),X_2(2))- eval_hat_functions(X_3(1),X_3(2),i,msh)).^2;
val=abs(detJ)*S*(1/3)*val;
end


