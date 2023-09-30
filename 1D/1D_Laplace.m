%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Teoría y Práctica Elementos Finitos
% Presenta: Luis Edwin Aguilar Anzures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descripcion: Dada f calculamos el error en L^2 y H^1 y sus respectivas
% tasas de decrecimiento para divisiones de 10, 20, 40, 80 y 160 entre ella
% misma y la solución aproximada obtenida con el elemento finito de lagrange de orden 1.
%Ouput: tabla de erroren L^2 y H^1 y sus tasas de decrecimiento 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%El programa hace una iteración para cada una de la divisiones de
%intervalos:10, 20, 40, 80 y 160


%% Parametros de Usuario
%interval endpoints of [a,b]
a=0;
b=1;
nI_approx_init=10;      % numero inicial del total de sub-intervalos en el mallado
n_cicles=5;             % numero de ciclos que corre el algoritmo numérico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Para las iteraciones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
nI_approx_vec=zeros(1,n_cicles); %Construimos un vector que obtiene el numero total de puntos en el mallado
  %para cada ciclo de refinamiento 
nI_approx_vec(1)=nI_approx_init; %El primer ciclo no se refina
for k=2:n_cicles
  nI_approx_vec(k)=2*nI_approx_vec(k-1); % El siguiente contiene el doble de puntos que el anterior
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2_error_vec=zeros(1,n_cicles); % Vector que contiene el error en L2 para cada ciclo
H1_error_vec=zeros(1,n_cicles); % Vector que contiene el error en H1 para cada ciclo

for k=1:n_cicles %aquí inicia el ciclo de iteraciones para cada refinamiento.

nI_approx=nI_approx_vec(k);        % numero de subintervalos  en el mallado para aproximar
nodes=linspace(a,b,nI_approx+1); %nodes de la malla, vector fila
n_dofs=length(nodes); % numero de DOFs (degrees of freedom)
L2_error_f=@L2_error;  %función que calcula el error L2 (ver descripción de la función para mas detalle)
H1_error_f=@H1_error; %función que calcula el error H1 (ver descripción de la función para mas detalle)


%para aproximar a u el subespacio finito resolvemos el sistema A u= b donde
%A es la matriz de difusión, b depende la los datos iniciales f y u son las coordenadas de la aproximación en el 
%espacio finito 

%%%%%%%%%%%%%%%%%%% CONTRUIMOS EL VECTOR b %%%%%%%%%%%%%%%%%%%
b_global=zeros(nI_approx+1,1); %es el vector b
b_local=zeros(2,1); %almacena las contribuciones de las integrales en los intervalos

  phi_m_f=@phi_m; %calcula la integral de f phi^-_i en (1,i+1)
  phi_p_f=@phi_p; %calcula la integral de f phi^+_{i+1} en (1,i+1)
for i=1:nI_approx
    b_local(1)= phi_m_f(nodes,i);
    b_local(2)= phi_p_f(nodes,i);
    b_global(i)= b_global(i)+b_local(1); %suma las respectivas contribuciones
    b_global(i+1)= b_global(i+1)+b_local(2); %suma las respectivas contribuciones
end
b_global(1)=0; %condiciones de frontera en a
b_global(nI_approx+1)=0; %condiciones de frontera en b


%%%%%%%%%%%%%%%%%%% CONTRUIMOS LA MATRIZ A  %%%%%%%%%%%%%%%%%%%
A_global=zeros(nI_approx+1,nI_approx+1); %es la matriz A
A_local=zeros(2,2);

 phi_nn_f=@phi_nn; %calcula la integral de (phi^-_i)'(phi^-_i)' en (1,i+1)
 phi_np_f=@phi_np;%calcula la integral de  (phi^-_i)'(phi^+_{i+1})' en (1,i+1)
 phi_pp_f=@phi_pp;%calcula la integral de  en (phi^+_{i+1})'(phi^+_{i+1})' en (1,i+1)
for i=1:nI_approx
     A_local(1,1)=phi_nn_f(nodes,i); %almacenamos las contribuciuones de las integrales
     A_local(2,2)=phi_pp_f(nodes,i);
     A_local(1,2)=phi_np_f(nodes,i);
     A_local(2,1)=phi_np_f(nodes,i);
%sumanos las contribuciones
    A_global(i,i)= A_global(i,i)+ A_local(1,1); 
    A_global(i+1,i+1)= A_global(i+1,i+1)+A_local(2,2);
    A_global(i,i+1)=A_local(1,2); 
    A_global(i+1,i)=A_local(2,1);
end
A_global(1,:)=0; %condicones de frontera en a
A_global(1,1)=1;
A_global(nI_approx+1,:)=0; %condicones de frontera en b
A_global(nI_approx+1,nI_approx+1)=1;

%%%%%%%%%%%%%%%%%%% RESOLVEMOS PARA u %%%%%%%%%%%%%%%%%%%%%
dofs=linsolve(A_global,b_global); %obtenemos el vector u de grados de libertad solución del sistema A|u_h=b
%los grados de libertad  son las coordenadas de la aproximación en el 
%espacio finito generado por las hat functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculo del eror u - u_h en la norma L2
    error=0;
    %Ciclo para cada subintervalo
    for i=1:nI_approx  
        error= error+ L2_error_f(nodes, dofs,i); %suma el valor de la integral en cada subintervalo 
    end
    error= sqrt(error);
    L2_error_vec(k)=error; %Guardar el valor del error L2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculo del eror u - u_h en la norma H1
   h1error=0;
    %Ciclo para cada subintervalo
    for i=1:nI_approx  
        h1error= h1error+ H1_error_f(nodes, dofs,i);  %suma el valor de la integral en cada subintervalo 
    end
    h1error= sqrt(h1error);
    H1_error_vec(k)=h1error; %Guardar el valor del error H1
    
end %acaba el ciclo de refinamientos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Crear tabla con valores: N| error L2| tasa_de_decrecimiento L2| error H1|
% tasa_de_decrecimiento H1
disp("##### Final del algoritmo #########"); %Imprime en la terminal 
% decrecimiento del error L2
err_rate=zeros(1,n_cicles); %tasa de decrecimiento del error
err_rate(1)=1;
for i=2:n_cicles
    err_rate(i)=log(L2_error_vec(i)./L2_error_vec(i-1))./log(1/2);
end
% decrecimiento del error H1
H1err_rate=zeros(1,n_cicles); %tasa de decrecimiento del error
H1err_rate(1)=1;
for i=2:n_cicles
    H1err_rate(i)=log(H1_error_vec(i)./H1_error_vec(i-1))./log(1/2);
end


%renombramos lo títulos de la tabla: 
N_vec=nI_approx_vec';
error_L2=L2_error_vec';
radio_de_error_L2= err_rate';
error_H1=H1_error_vec';
radio_de_error_H1= H1err_rate';


output_table=table(N_vec,error_L2,radio_de_error_L2,error_H1,radio_de_error_H1) ; % Crea tabla con vectores columna                                           
disp(output_table);%Imprime en la terminal

%%%%%%%%%%%%%% FUNCIONES %%%%%%%%%%%%%%%%%

%la sol teórica 
function [val] = u_exact(x)
    val = sin(4*pi*x); 
end
%la derivada de la sol teórica
function [val] = du_exact(x)
    val = 4*pi*cos(4*pi*x); 
end

% Function:phi_m
%       nodes:  los nodos de la particion [a,b]
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       Evalua la integral (\phi^-_i)f  en el intervalo  (nodes(i),nodes(i+1))
%      usando la regla de(l) Simpson 
function [val] = phi_m(nodes,i)
            x0 =nodes(i);
            x2 =nodes(i+1);
             h=0.5*(x2-x0);
            x1 = x0+h;
            val= ((pi*4)^2)*sin(4*pi*x0)*(x0-x2)/(x0-x2);
            val= val +4* ((pi*4)^2)*sin(4*pi*x1)*(x1-x2)/(x0-x2);
            val= val + ((pi*4)^2)*sin(4*pi*x2)*(x2-x2)/(x0-x2);
            val = h/3.*val;
end

% Function:phi_p
%       nodes:  los nodos de la particion [a,b]
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       Evalua la integral (\phi^+_{i+1})f  en el intervalo  (nodes(i),nodes(i+1))
%       usando la regla de(l) Simpson
function [val] = phi_p(nodes,i)
            x0 =nodes(i);
            x2 =nodes(i+1);
            h=0.5*(x2-x0);
            x1 = x0+h;
            val= ((pi*4)^2)*sin(4*pi*x0)*(x0-x0)/(x2-x0);
            val= val +4* ((pi*4)^2)*sin(4*pi*x1)*(x1-x0)/(x2-x0);
            val= val + ((pi*4)^2)*sin(4*pi*x2)*(x2-x0)/(x2-x0);
            val = h/3.*val;
end

%%%% el siguiente grupo de funciones es para calcular las integrales en la
%%%% matriz de difusión, en este caso aunque salen constantes aplicamos la
%%%% regla de Simpson en aras de la numericidad 


% Function:phi_nn
%       nodes:  los nodos de la particion [a,b]
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       alcula la integral de (phi^-_i)'(phi^-_i)' en el intervalo (1,i+1)
%      usando la regla de(l) Simpson 
function [val] = phi_nn(nodes,i)
            x0 =nodes(i);
            x2 =nodes(i+1);
             h=0.5*(x2-x0);
            val= 1/((x0-x2)^2);
             val=val+ 4/((x0-x2)^2);
             val=val+ 1/((x0-x2)^2);
            val = h/3.*val;
end

% Function:phi_np
%       nodes:  los nodos de la particion [a,b]
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       alcula la integral de  (phi^-_i)'(phi^+_{i+1})'  en el intervalo (1,i+1)
%      usando la regla de(l) Simpson 
function [val] = phi_np(nodes,i)
            x0 =nodes(i);
            x2 =nodes(i+1);
             h=0.5*(x2-x0);
            val= (1/(x0-x2)) * (1/(x2-x0));
            val=val + 4*(1/(x0-x2))*(1/(x2-x0));
            val=val + (1/(x0-x2))*(1/(x2-x0));
            val = h/3.*val;
end


% Function:phi_pp
%       nodes:  los nodos de la particion [a,b]
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       alcula la integral de (phi^+_{i+1})'(phi^+_{i+1})'  en el intervalo (1,i+1)
%      usando la regla de(l) Simpson 
function [val] = phi_pp(nodes,i)
            x0 =nodes(i);
            x2 =nodes(i+1);
             h=0.5*(x2-x0);
            val= 1/((x2-x0)^2);
             val=val+ 4/((x2-x0)^2);
             val=val+ 1/((x2-x0)^2);
            val = h/3.*val;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Functio: eval_hat_function_pair
% Input: 
%       x:      punto a evaluar
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes
%       i      : el índice de la hat_función  de la izquierda
%              
% Descripcion:
%               Evalua en el punto x el par de funciones con índices 
%               i e i+1 utilizando a su vez los  los grados de libertad
%               dofs.
%               Precaucion: Para hacer esta función mucho más compacta
%               se asume que x efectivamente esta en el
%               intervalo [nodes(i), nodes(i+1)]
function [val] = eval_hat_function_pair(x,nodes, dofs,i)

        eps=1e-6; %Tolerancia
        %Checamos si x esta en el intervalo donde hat_i(x) es cero
        if((x+eps<nodes(i)) || (x-eps > nodes(i+1))) % Ver Gilat sec. 6.1 & 6.2
           disp("Error: x is not defined in hat_functions_pair ") ;%Gilat sec 4.3.1
           assert(false); %Stops program
        end
       val= dofs(i)*(1- (x- nodes(i))/(nodes(i+1)-nodes(i)));
       val= val + dofs(i+1)*(x-nodes(i))/(nodes(i+1)-nodes(i));
       
    
end


% Functio: eval_diff_hat_function_pair
% Input: 
%       x:      punto a evaluar
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes
%       i      : el índice de la hat_función  de la izquierda
%              
% Descripcion:
%               Evalua en el punto x eln la derivada del par de funciones con índices 
%               i e i+1 utilizando a su vez los  los grados de libertad
%               dofs.
%               Precaucion: Para hacer esta función mucho más compacta
%               se asume que x efectivamente esta en el
%               intervalo [nodes(i), nodes(i+1)]
function [val] = eval_diff_hat_function_pair(x,nodes, dofs,i)

        eps=1e-6; %Tolerancia
        %Checamos si x esta en el intervalo donde hat_i(x) es cero
        if((x+eps<nodes(i)) || (x-eps > nodes(i+1))) % Ver Gilat sec. 6.1 & 6.2
           disp("Error: x is not defined in hat_functions_pair ") ;%Gilat sec 4.3.1
           assert(false); %Stops program
        end
 
      val= dofs(i)*(1/(nodes(i)-nodes(i+1)));
       val= val + dofs(i+1)*(1/(nodes(i+1)-nodes(i)));
    
end


% Function: L2-error
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes       
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       Evalua la integral (u-uh)^2 en el intervalo  (nodes(i),nodes(i+1))
%       utilizando  la regla del Simpson. Ver Libro de Burden-Faires pg.196 
function [val] = L2_error(nodes, dofs,i)
            x0 =nodes(i);
            x2 =nodes(i+1);
            h=0.5*(x2-x0);
            x1 = x0+h;
          
            val= (u_exact(x0) - eval_hat_function_pair(x0,nodes, dofs,i))^2;
            val= val + 4*(u_exact(x1) - eval_hat_function_pair(x1,nodes, dofs,i))^2;
            val= val + (u_exact(x2) - eval_hat_function_pair(x2,nodes, dofs,i))^2;
            val = h/3.*val;
end

% Function: sem_H1_error
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes       
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       Evalua la integral (u'-uh')^2 en el intervalo  (nodes(i),nodes(i+1))
%       utilizando  la regla del Simpson. Ver Libro de Burden-Faires pg.196 
%la integral en cuentión corresponde a la seminorma de Sobolev en dicho
%intervalo

function [val] = sem_H1_error(nodes, dofs,i)
            x0 =nodes(i);
            x2 =nodes(i+1);
            h=0.5*(x2-x0);
            x1 = x0+h;
          
            val= (du_exact(x0) - eval_diff_hat_function_pair(x0,nodes, dofs,i))^2;
            val= val + 4*(du_exact(x1) - eval_diff_hat_function_pair(x1,nodes, dofs,i))^2;
            val= val + (du_exact(x2) - eval_diff_hat_function_pair(x2,nodes, dofs,i))^2;
            val = h/3.*val;
end

% Function: H1_error
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes       
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       Evalua la integral (u'-uh')^2 + (u-uh)^2 en el intervalo  (nodes(i),nodes(i+1))
%       utilizando  las funciones sem_H1_error(nodes, dofs,i) y L2_error(nodes, dofs,i)

function [val] = H1_error(nodes, dofs,i)
          val= L2_error(nodes, dofs,i)+sem_H1_error(nodes, dofs,i);
end


