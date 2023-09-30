%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Teoría y Práctica Elementos Finitos
% Presenta: Luis Edwin Aguilar Anzures
% email: eda.math10@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descripcion: Aproximámos numéricamente utilizando los elementos finitos de Lagrange de primer orden 
% la solució débil del problema dado.
% Ouput: Obtenemos la tasa de convergencia para el error en norma L^2 y
% norma H^1 entre la solución teórica y la solución aproximada.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cargamos las mallas dadas
msh_filenames = ["square_1.m" "square_2.m" "square_3.m"...
  "square_4.m" "square_5.m" "square_6.m"];
h_vec=zeros(1,6);% Vector que contiene el h para cada malla
L2_error_vec=zeros(1,6); % Vector que contiene el error en L2 para cada malla
H1_error_vec=zeros(1,6); % Vector que contiene el error en H1 para cada ciclo
u_x_y=@u_exact; %cargamos la solución teórica u


%empieza el ciclo para cada malla dada.
for j=1:6
msh=read_mesh(msh_filenames(j));
h_vec(j)=get_h_global(msh);
%%%%%%%%%%%%%%%%%%% CONSTRUCCIÓN DEL SISTEMA DE GALERKIN %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%

a_bilear= @a_bilienal_global;%función que calcula las contribuciones de la forma bilineal de las hat functions globales en un elemento. Ver descripción de la función para más infiormación. 
int_f_phi =@int_fphi; %función que calcula las contribuciones de ña integral de f con las hat functions globales en un elelemnnto. Ver descripción de la función para más infiormación. 
A_g=zeros(msh.nb_nodes,msh.nb_nodes); %donde se guardaran los datos de la matriz asociada al sistema
b_g=zeros(msh.nb_nodes,1); %donde se guardaran los datos del vector asociado asociado a la función f del sistema

%%%%%%%%%%%%%%%%%%% ENSAMBLAJE DEL SISTEMA%%%%%%%%%%%%%%%%%%%
% Siguiendo el código de Finite Elements and Fast Iterative Solvers with Applications in 
% Incompressible Fluid Dynamics, H. Elman, D. Silvester, and A. Wathen.,
% (2005) en la parte 1.4.3 Assembly of the Galerkin system
for s = 1:msh.nb_elems %empieza el ciclo por elemento
   for k = 1:3
  for i = 1:3%se suman las conbtribuciones de las iontegrales asociadas a cada nodo que yace en ele elemento fijado 
A_g( msh.elems_nodes_conn(s,i),msh.elems_nodes_conn(s,k)) = A_g( msh.elems_nodes_conn(s,i),msh.elems_nodes_conn(s,k)) + a_bilear(msh.elems_nodes_conn(s,i),msh.elems_nodes_conn(s,k),s,msh);
  end %termina ciclo pra i
    b_g(msh.elems_nodes_conn(s,k)) = b_g(msh.elems_nodes_conn(s,k)) +  int_f_phi(msh.elems_nodes_conn(s,k),s,msh) ;
   end %termina ciclo pra k
end %termina ciclo pra s

%Condiciones de frontera: En la numeración de los nodos, todos los nodos de
%la frontera se encuentran al inicio y tomamos ventaja de este hecho;
 indicesCeros = find(msh.onodes_bool== 0); %hallamos el primer nodo que no sea de forntera 
%a partir de este indice todos son nodos interiores
     primer_int = indicesCeros(1);
for i=1: primer_int-1
   %las siguientes dos condiciones capturan el hecho de que los nodos de
   %fronetra aportan a la integral total.
A_g(i,:)=zeros(msh.nb_nodes,1)  ;
A_g(:,i)=zeros(1,msh.nb_nodes) ;
 %las siguientes dos condiciones son para restringirnos al espacio de hat
 %function globales qyue valen cero en la frontera, es decir para que se
 %anulen las hat functión asociadas a nodos de frontera (pues asi es el
 %espacio finito en el que aproximamos).
A_g(i,i)=1;
b_g(i)=0;
end    
 


%%%%%%%%%%%%%%%%%%% RESOLVEMOS PARA U %%%%%%%%%%%%%%%%%%%%%
dofs=linsolve(A_g,b_g); %obtenemos el vector u de grados de libertad solución del sistema A|u_h=b
%los grados de libertad  son las coordenadas de la aproximación en el 
%espacio finito generado por las hat functions

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculo del eror u - u_h en la norma L2

%sumamos los errores L^2 en cada elelemnto para obtener el error total de
%la malla 

    L2_error_f= @L2_error; %función que calcula el error en cada elemento. Vea descripción de la función.
    error=0;
    %Ciclo para cad amalla
    for i=1:msh.nb_elems  
       error=error+ L2_error_f(i,dofs,msh);
    end
    error= sqrt(error);
    L2_error_vec(j)=error; %Guardar el valor del error L2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Calculo del eror u - u_h en la norma H1

%sumamos los errores de la seminorma de Sobolev H^1 en cada elemento para obtener el error total de
%la malla y al final le sumamos el respectivo error en L^2 para calcular el
%error en la norma de Sobolev H^1.

H1_error_f= @semH1_error ; %función que calcula el error en seminorma en cada elemento. Vea descripción de la función.
   h1error=0;
    %Ciclo para cada subintervalo
    for i=1:msh.nb_elems  
        h1error= h1error+ H1_error_f(i,dofs,msh);  %suma el valor de la integral en cada subintervalo 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%la sol teórica 
function [val] = u_exact(x,y)
    val = sin(x*pi*4).*sin(y*pi*4); 
end
%la derivada parcial respecto a x sol teórica
function [val] = dx_u_exact(x,y)
    val = pi*4*cos(x*pi*4).*sin(y*pi*4); 
end
%la derivada parcial respecto a y sol teórica
function [val] = dy_u_exact(x,y)
    val =pi*4*sin(x*pi*4).*cos(y*pi*4); 
end

%función de la forma bilinal
function [val] = K(x)
    val = 1+x^2;
end
%función fuente 
function [val] = f(x,y)
    val = 8*pi*sin(4*pi*y)*(4*pi*(x^2)*sin(4*pi*x)+4*pi*sin(4*pi*x)-x*cos(4*pi*x));
end


%%DESCRIPCIÓN: Calcularmos el interpolador de en un elemento i usando las funciones de
%%forma del elemento de referencia y una transformación al elemento de
%%referencia 
function [val] = eval_hat_functions(x,y,dofs,i,msh)
%msh=read_mesh(['square_',num2str(j),'.m']);
V_1=msh.nodes(msh.elems_nodes_conn(i,1),1:2); %primer vértice del elemnto i
V_2=msh.nodes(msh.elems_nodes_conn(i,2),1:2); %segundo vértice del elemnto i
V_3=msh.nodes(msh.elems_nodes_conn(i,3),1:2); %tercer vértice del elemnto i
IN=inpolygon(x,y,[V_1(1),V_2(1),V_3(1)],[V_1(2),V_2(2),V_3(2)]);  %revisa si el punto (x,y) yace el elemento i
A=[V_3'-V_1',V_2'-V_1']; %matriz del mapeo del elemento de refenrencia al elemento i  T(x,y)=A(x,y)+V_1
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)]; %matriz del mapeo del elemento el elemento i al de referencia  T(x,y)=A_inv((x,y)-V_1)
val=dofs(msh.elems_nodes_conn(i,1))*(1-(x-V_1(1))*A_in(1,1)-(y-V_1(2))*A_in(1,2)-(x-V_1(1))*A_in(2,1)-(y-V_1(2))*A_in(2,2) ); %evalucación en la función de forma del vértice V_1
val=val+ dofs(msh.elems_nodes_conn(i,2))*((x-V_1(1))*A_in(1,1)+(y-V_1(2))*A_in(1,2)); %suma evalucación en la función de forma del vértice V_2
val=val+ dofs(msh.elems_nodes_conn(i,3))*((x-V_1(1))*A_in(2,1)+(y-V_1(2))*A_in(2,2)); %evalucación en la función de forma del vértice V_3
val=IN.*val; %asegura que el soporte de este interpolador se quede en el elemnto i
end



%%DESCRIPCIÓN: Calcularmos la derivada parcial respecto a la coordenada x interpolador(x,y) de en un elemento i usando las funciones de
%%forma del elemento de referencia y una transformación al elemento de
%%referencia 
function [val] = dx_eval_hat_functions(x,y,dofs,i,msh)
V_1=msh.nodes(msh.elems_nodes_conn(i,1),1:2);
V_2=msh.nodes(msh.elems_nodes_conn(i,2),1:2);
V_3=msh.nodes(msh.elems_nodes_conn(i,3),1:2);
IN=inpolygon(x,y,[V_1(1),V_2(1),V_3(1)],[V_1(2),V_2(2),V_3(2)]) ;
A=[V_3'-V_1',V_2'-V_1'];
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)];
val=dofs(msh.elems_nodes_conn(i,1))*(-A_in(1,1)-A_in(2,1));
val = val +dofs(msh.elems_nodes_conn(i,2))*A_in(1,1);
val = val +dofs(msh.elems_nodes_conn(i,3))*A_in(2,1);
val=IN.*val;
end

%%DESCRIPCIÓN: Calcularmos la derivada parcial respecto a la coordenada y interpolador(x,y) de en un elemento i usando las funciones de
%%forma del elemento de referencia y una transformación al elemento de
%%referencia 
function [val] = dy_eval_hat_functions(x,y,dofs,i,msh)
V_1=msh.nodes(msh.elems_nodes_conn(i,1),1:2);
V_2=msh.nodes(msh.elems_nodes_conn(i,2),1:2);
V_3=msh.nodes(msh.elems_nodes_conn(i,3),1:2);
IN=inpolygon(x,y,[V_1(1),V_2(1),V_3(1)],[V_1(2),V_2(2),V_3(2)]) ;
A=[V_3'-V_1',V_2'-V_1'];
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)];
val=dofs(msh.elems_nodes_conn(i,1))*(-A_in(1,2)-A_in(2,2));
val = val +dofs(msh.elems_nodes_conn(i,2))*A_in(1,2);
val = val +dofs(msh.elems_nodes_conn(i,3))*A_in(2,2);
val=IN.*val;
end

%DESCRIPCIÓN: función que mapea un punto en el elementos de referencia al elemento i
%para hacer la integración T(x,y)=A(x,y)+V_1, donde A=[V_3-_V_1,V_2-_V_1]
function [val] =map(x,y,i,msh)
%msh=read_mesh(['square_',num2str(j),'.m']);
val=[x*(msh.nodes(msh.elems_nodes_conn(i,3),1)-msh.nodes(msh.elems_nodes_conn(i,1),1))+y*(msh.nodes(msh.elems_nodes_conn(i,2),1)-msh.nodes(msh.elems_nodes_conn(i,1),1))+msh.nodes(msh.elems_nodes_conn(i,1),1),x*(msh.nodes(msh.elems_nodes_conn(i,3),2)-msh.nodes(msh.elems_nodes_conn(i,1),2))+y*(msh.nodes(msh.elems_nodes_conn(i,2),2)-msh.nodes(msh.elems_nodes_conn(i,1),2))+msh.nodes(msh.elems_nodes_conn(i,1),2) ];
end


 %DESCRIPCIÓN:Calculla el error en la seminorma de Sobolvev en el elemnto i usando un cambio de variable al elemento 
 % de referencia y el método de cuadratura Gaussiana de orden dos (tres puntos)
function [val]=semH1_error(i,dofs,msh)
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
val= (dx_u_exact(X_1(1),X_1(2))- dx_eval_hat_functions(X_1(1),X_1(2),dofs,i,msh)).^2+(dy_u_exact(X_1(1),X_1(2))- dy_eval_hat_functions(X_1(1),X_1(2),dofs,i,msh)).^2;
  val= val + (dx_u_exact(X_2(1),X_2(2))- dx_eval_hat_functions(X_2(1),X_2(2),dofs,i,msh)).^2+(dy_u_exact(X_2(1),X_2(2))- dy_eval_hat_functions(X_2(1),X_2(2),dofs,i,msh)).^2;
  val= val + (dx_u_exact(X_3(1),X_2(2))- dx_eval_hat_functions(X_3(1),X_3(2),dofs,i,msh)).^2+(dy_u_exact(X_3(1),X_2(2))- dy_eval_hat_functions(X_3(1),X_3(2),dofs,i,msh)).^2;
val=abs(detJ)*S*(1/3)*val;
end



%DESCRIPCIÓN:Calculla el error en la seminorma de Sobolvev en el elemnto i usando un cambio de variable al elemento 
 % de referencia y el método de cuadratura Gaussiana de orden dos (tres puntos)
function [val] =L2_error(i,dofs,msh)

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
val= (u_exact(X_1(1),X_1(2))- eval_hat_functions(X_1(1),X_1(2),dofs,i,msh)).^2;
  val= val + (u_exact(X_2(1),X_2(2))- eval_hat_functions(X_2(1),X_2(2),dofs,i,msh)).^2;
  val= val + (u_exact(X_3(1),X_2(2))- eval_hat_functions(X_3(1),X_3(2),dofs,i,msh)).^2;
val=abs(detJ)*S*(1/3)*val;
end



%NOTA: Cuando digo indices locales me refiero a los indices de los nodos
%relativos al elemento donde se trabaja. Es decir, Fijo un elemento con
%índice s, los índies de  los nodos locales son: 1, 2 y 3. Esto lo hago
%así, pues para ensamblar el sistema de Galerkin vamos sumando las contribuciones por de
%las inbtegrales por cada elemento 



%DESCRIPCIÓN: Aplica la forma bileneal de la formulación débil a las hat
%functions globales con índices locales i y k (centradas en i y k respectivamente) en el elemento con índice s
%i.e., calcula la contribución de a(\phi_i,\phi_k) en el elemento con
%índice s. Esto lo hace cuando i es distinto de k.
 function val = a_bilear_int(i,k,s,msh)
otro_vertice=setdiff( msh.elems_nodes_conn(s,1:3), [i,k] ) ;
V_3=msh.nodes(i,1:2); %nodo local i 
V_2=msh.nodes(k,1:2); %nodo local k 
V_1=msh.nodes(otro_vertice,1:2); %el nodo restante.
A=[V_3'-V_1',V_2'-V_1']; %matriz asociada al mapeo del elemento de referencia al elemento con índice s
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)]; %matriz asociada al mapeo del elemento con índice sal al elemento de referencia 
dx_lamnda_2 = A_in(1,1); %parcial respecto a x de la hat function asociada al nodo local i 
dx_lamnda_3 = A_in(2,1); %parcial respecto a x de la hat function asociada al nodo local k 
dy_lamnda_2 = A_in(1,2); %parcial respecto a y de la hat function asociada al nodo local i 
dy_lamnda_3 = A_in(2,2);%parcial respecto a y de la hat function asociada al nodo local k
S =msh.nodes(msh.elems_nodes_conn(s,1),1)*(msh.nodes(msh.elems_nodes_conn(s,2),2)-msh.nodes(msh.elems_nodes_conn(s,3),2))+msh.nodes(msh.elems_nodes_conn(s,2),1)*(msh.nodes(msh.elems_nodes_conn(s,3),2)-msh.nodes(msh.elems_nodes_conn(s,1),2))+msh.nodes(msh.elems_nodes_conn(s,3),1)*(msh.nodes(msh.elems_nodes_conn(s,1),2)-msh.nodes(msh.elems_nodes_conn(s,2),2));
 S=0.5*S; %área del elemnto s 
X_1=(A*[0,0.5]'+V_1')'; %primer punto para la cuadratura
X_2=(A*[0.5,0.5]'+V_1')';  %segundo punto para la cuadratura
X_3=(A*[0.5,0]'+V_1')';%tercer punto para la cuadratura
% X_1=map(0,0.5,s,msh); %primer punto para la cuadratura
% X_2=map(0.5,0.5,s,msh);  %segundo punto para la cuadratura
% X_3=map(0.5,0,s,msh);%tercer punto para la cuadratura
%determinante de la derivada del mapeo que va del simplejo de referencia al
%elekemnto s 
detJ=det([msh.nodes(msh.elems_nodes_conn(s,3),1:2)'-msh.nodes(msh.elems_nodes_conn(s,1),1:2)',msh.nodes(msh.elems_nodes_conn(s,2),1:2)'-msh.nodes(msh.elems_nodes_conn(s,1),1:2)'] ) ;
%evaluamos la integral de acuerdo con la fórmula de la cuadratura
%gaussiana a tres puntos.
val= K(X_1(1))*((dx_lamnda_3*dx_lamnda_2)+(dy_lamnda_3*dy_lamnda_2)) ;
val=val + K(X_2(1))*((dx_lamnda_3*dx_lamnda_2)+(dy_lamnda_3*dy_lamnda_2)) ;
val=val + K(X_3(1))*((dx_lamnda_3*dx_lamnda_2)+(dy_lamnda_3*dy_lamnda_2)) ;
val=abs(detJ)*S*(1/3)*val;
end

%DESCRIPCIÓN: Aplica la forma bileneal de la formulación débil a hat
%function global con índices locales i (centrada en i) con ella misma en el elemento con índice s
%i.e., calcula la contribución de a(\phi_i,\phi_i) en el elemento con
%índice s. Estos valores corresponadran a la diagonal de la matriz del sistema de Galerkin 
 function val = a_bilear_int_diag(i,s,msh)
otro_vertice=setdiff( msh.elems_nodes_conn(s,1:3), i ) ;
V_3=msh.nodes(i,1:2); %nodo local i 
V_2=msh.nodes(otro_vertice(2),1:2); %los nodos que restan
V_1=msh.nodes(otro_vertice(1),1:2); 
A=[V_3'-V_1',V_2'-V_1'];%matriz asociada al mapeo del elemento de referencia al elemento con índice s
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)];%matriz asociada al mapeo del elemento con índice sal al elemento de referencia 
dx_lamnda_3 = A_in(2,1);%parcial respecto a x de la hat function asociada al nodo local i 
dy_lamnda_3 = A_in(2,2);%parcial respecto a y de la hat function asociada al nodo local i 
S =msh.nodes(msh.elems_nodes_conn(s,1),1)*(msh.nodes(msh.elems_nodes_conn(s,2),2)-msh.nodes(msh.elems_nodes_conn(s,3),2))+msh.nodes(msh.elems_nodes_conn(s,2),1)*(msh.nodes(msh.elems_nodes_conn(s,3),2)-msh.nodes(msh.elems_nodes_conn(s,1),2))+msh.nodes(msh.elems_nodes_conn(s,3),1)*(msh.nodes(msh.elems_nodes_conn(s,1),2)-msh.nodes(msh.elems_nodes_conn(s,2),2));
 S=0.5*S; %área del elemnto s 
 X_1=map(0,0.5,s,msh); %primer punto para la cuadratura
 X_2=map(0.5,0.5,s,msh);  %segundo punto para la cuadratura
 X_3=map(0.5,0,s,msh);%tercer punto para la cuadratura
%determinante de la derivada del mapeo que va del simplejo de referencia al
%elekemnto s 
detJ=det([msh.nodes(msh.elems_nodes_conn(s,3),1:2)'-msh.nodes(msh.elems_nodes_conn(s,1),1:2)',msh.nodes(msh.elems_nodes_conn(s,2),1:2)'-msh.nodes(msh.elems_nodes_conn(s,1),1:2)'] ) ;
%evaluamos la integral de acuerdo con la fórmula de la cuadratura
%gaussiana
val= K(X_1(1))*((dx_lamnda_3)^2+(dy_lamnda_3)^2) ;
val=val + K(X_2(1))*((dx_lamnda_3)^2+(dy_lamnda_3)^2) ;
val=val + K(X_3(1))*((dx_lamnda_3)^2+(dy_lamnda_3)^2) ;
val=abs(detJ)*S*(1/3)*val;
end

%DESCRIPCIÓN:Calcula la contribución de la integral de f y la hat functión
%con índice local i (centrada en i) (relativa al elemento s) en el elemento
% con índice s
function val = int_fphi(i,s,msh)
otro_vertice=setdiff( msh.elems_nodes_conn(s,1:3), i ) ;
V_3=msh.nodes(i,1:2); %nodo local i 
V_2=msh.nodes(otro_vertice(2),1:2); %los otros nodos del elemento s
V_1=msh.nodes(otro_vertice(1),1:2);
A=[V_3'-V_1',V_2'-V_1'];%matriz asociada al mapeo del elemento de referencia al elemento con índice s
A_in=(1/det(A))*[V_2(2)-V_1(2),V_1(1)-V_2(1);V_1(2)-V_3(2),V_3(1)-V_1(1)];%matriz asociada al mapeo del elemento con índice sal al elemento de referencia 
S =msh.nodes(msh.elems_nodes_conn(s,1),1)*(msh.nodes(msh.elems_nodes_conn(s,2),2)-msh.nodes(msh.elems_nodes_conn(s,3),2))+msh.nodes(msh.elems_nodes_conn(s,2),1)*(msh.nodes(msh.elems_nodes_conn(s,3),2)-msh.nodes(msh.elems_nodes_conn(s,1),2))+msh.nodes(msh.elems_nodes_conn(s,3),1)*(msh.nodes(msh.elems_nodes_conn(s,1),2)-msh.nodes(msh.elems_nodes_conn(s,2),2));
 S=0.5*S; %área del elemnto s 
 X_1=map(0,0.5,s,msh); %primer punto para la cuadratura
 X_2=map(0.5,0.5,s,msh);  %segundo punto para la cuadratura
 X_3=map(0.5,0,s,msh);%tercer punto para la cuadratura
%determinante de la derivada del mapeo que va del simplejo de referencia al
%elekemnto s 
detJ=det([msh.nodes(msh.elems_nodes_conn(s,3),1:2)'-msh.nodes(msh.elems_nodes_conn(s,1),1:2)',msh.nodes(msh.elems_nodes_conn(s,2),1:2)'-msh.nodes(msh.elems_nodes_conn(s,1),1:2)'] ) ;
%evaluamos la integral de acuerdo con la fórmula de la cuadratura
%gaussiana
val= f(X_1(1),X_1(2))*((X_1(1)-V_1(1))*A_in(2,1)+(X_1(2)-V_1(2))*A_in(2,2));
val=val + f(X_2(1),X_2(2))*((X_2(1)-V_1(1))*A_in(2,1)+(X_2(2)-V_1(2))*A_in(2,2));
val=val + f(X_3(1),X_3(2))*((X_3(1)-V_1(1))*A_in(2,1)+(X_3(2)-V_1(2))*A_in(2,2)) ;
val=abs(detJ)*S*(1/3)*val;
end

%DESCRIPCIÓN: Aplica la forma bileneal de la formulación débil a las hat
%functions globales con índices locales i y k (centradas en i y k respectivamente) en el elemento con índice s
%i.e., calcula la contribución de a(\phi_i,\phi_k) en el elemento con
%índice s. Esto la hace para cualquier par de íncides locales i  y k (relativos a s) sin importar si son
%iguales 
function val=a_bilienal_global(i,k,s,msh)
  if i==k
   val=a_bilear_int_diag(i,s,msh); %si los indices locales son distintos 
  else 
   val=a_bilear_int(i,k,s,msh); %si los indices locales son iguales
  end
end



