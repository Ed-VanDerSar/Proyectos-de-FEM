%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Teoría y Práctica Elementos Finitos
% Posgrado en Matematicas-UNAM-CdMx
% Presenta: Luis Edwin Aguilar Anzures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descripcion:  Dada u calculamos el error en L^2 y H^1 entre u misma 
% y su interpolador respecto a los elementos finitos de lagrange de grado 2
% y sus respectivas tasas de decrecimiento para divisiones de 10, 20, 40, 80 y 160 subinterivalos
%Ouput: tabla de errores L^2 y H^1 y sus tasas de decrecimiento 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%El programa hace una iteración para cada una de la divisiones de
%intervalos:10, 20, 40, 80 y 160

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametros de Usuario
%interval endpoints of [a,b]
a=0;
b=1;
nI_approx_init=10;      % numero inicial del total de sub-intervalos en el mallado
n_cicles=5;    % numero de ciclos que corre el algoritmo numérico
u_exact_f=@u_exact; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Para las iteraciones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
nI_approx_vec=zeros(1,n_cicles); %Construimos un vector que obtiene el numero total de puntos en el mallado
  %para cada ciclo de refinamiento 
nI_approx_vec(1)=nI_approx_init; %El primer ciclo no se refina
for k=2:n_cicles
  nI_approx_vec(k)=2*nI_approx_vec(k-1); % El siguiente contiene el doble de puntos que el anterior
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2_error_vec=zeros(1,n_cicles); % Vector que contiene el error en L2 para cada ciclo
H1_error_vec=zeros(1,n_cicles); % Vector que contiene el error en L2 para cada ciclo

for k=1:n_cicles

nI_approx=nI_approx_vec(k);        % numero de subintervalos  en el mallado para aproximar
nodes=linspace(a,b,nI_approx+1); %nodes de la malla, vector fila
L2_error_f=@L2_error;  % Vector que contiene el error en L2 para cada ciclo
H1_error_f=@H1_error;  % Vector que contiene el error en H1 para cada ciclo
int_nodes=linspace(a,b,2*nI_approx+1); %nodos para calcular el error en las integrales 
dofs=u_exact_f(int_nodes)';  %vector columna (transpuesto)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Calculo del eror u - u_h en la norma L2
    error=0;
    %Ciclo para cada subintervalo
    for i=1:nI_approx  
        error= error+ L2_error_f(nodes, dofs,i);
    end
    error= sqrt(error);
    L2_error_vec(k)=error; %Guardar el valor del error L2


  %  nodes_int=linspace(a,b,2*nI_approx+1):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculo del eror u - u_h en la norma H1
   h1error=0;
    %Ciclo para cada subintervalo
    for i=1:nI_approx  
        h1error= h1error+ H1_error_f(nodes, dofs,i);
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
    val = cos(4*pi*x); 
end
%la derivada de la sol teórica
function [val] = du_exact(x)
    val = -4*pi*sin(4*pi*x); 
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
   r=(nodes(i+1)-nodes(i))^2 ;
    val=dofs(i+i-1)*(2*x^2 -x*(nodes(i)+3*nodes(i+1))+ nodes(i+1)*(nodes(i)+nodes(i+1))) ;
    val=val+dofs(i+i+1)*(2*x^2 -x*(nodes(i+1)+3*nodes(i))+ nodes(i)*(nodes(i)+nodes(i+1))) ;
    val=val+ dofs(i+i)*4*(x*(nodes(i+1)+nodes(i))- (nodes(i)*nodes(i+1))-x^2) ;
    val=(1/r) * val  ;
   % val=dofs(i+1)*(x-nodes(i))*(x-pm)/((nodes(i+1)-nodes(i))*(nodes(i+1)-pm)) ;
%val=val+dofs(i)*(x-nodes(i+1))*(x-pm)/((nodes(i)-nodes(i+1))*(nodes(i)-pm)) ;
%val=val+u_exact(pm)*(x-nodes(i+1))*(x-nodes(i))/((pm-nodes(i))*(pm-nodes(i+1)));
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
r=(nodes(i+1)-nodes(i))^2;
      val=dofs(i+i-1)*(4*x -(nodes(i)+3*nodes(i+1)));
      val=val+dofs(i+i+1)*(4*x -(nodes(i+1)+3*nodes(i)) );
      val=val+4*dofs(i+i)*(-2*x +(nodes(i+1)+nodes(i)));
    val=(1/r) * val  ;
end


% Function: L2-error
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes       
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       Evalua la integral (u-uh)^2 en el intervalo  (nodes(i),nodes(i+1))
%       utilizando  la regla del Simpson en cada subintervalo
%       ((nodes(i),pm)) y(pm,nodes(i+1)), donde pm es el punto medio 
function [val] = L2_error(nodes, dofs,i)
           h=(nodes(i+1)-nodes(i))/4;
           x0 =nodes(i);
            x1 = x0+h;
            x2 =(nodes(i+1)+nodes(i))/2;
            y0 =x2;
            y2 =nodes(i+1);
            y1 = (y0+y2)/2;
            
            val= (u_exact(x0) - eval_hat_function_pair(x0,nodes, dofs,i))^2;
            val= val + 4*(u_exact(x1) - eval_hat_function_pair(x1,nodes, dofs,i))^2;
            val= val + (u_exact(x2) - eval_hat_function_pair(x2,nodes, dofs,i))^2;
            val= val+(u_exact(y0) - eval_hat_function_pair(y0,nodes, dofs,i))^2;
            val= val + 4*(u_exact(y1) - eval_hat_function_pair(y1,nodes, dofs,i))^2;
            val= val + (u_exact(y2) - eval_hat_function_pair(y2,nodes, dofs,i))^2;
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
             h=(nodes(i+1)-nodes(i))/4;
           x0 =nodes(i);
            x1 = x0+h;
            x2 =(nodes(i+1)+nodes(i))/2;
            y0 =x2;
            y2 =nodes(i+1);
            y1 = (y0+y2)/2;
          
            val= (du_exact(x0) - eval_diff_hat_function_pair(x0,nodes, dofs,i))^2;
            val= val + 4*(du_exact(x1) - eval_diff_hat_function_pair(x1,nodes, dofs,i))^2;
            val= val + (du_exact(x2) - eval_diff_hat_function_pair(x2,nodes, dofs,i))^2;
            val= val+(du_exact(y0) - eval_diff_hat_function_pair(y0,nodes, dofs,i))^2;
            val= val + 4*(du_exact(y1) - eval_diff_hat_function_pair(y1,nodes, dofs,i))^2;
            val= val + (du_exact(y2) - eval_diff_hat_function_pair(y2,nodes, dofs,i))^2;
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
