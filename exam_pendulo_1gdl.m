%Limpieza de pantalla
clear all
close all
clc

%Declaración de variables simbólicas
     %theta1     %Velocidad angular     #Aceleración angular
syms th1(t)      th1p(t)                th1pp(t) 
syms m1 Ixx1 Iyy1 Izz1 %Masas y matrices de Inercia
syms l1 lc1   %l=longitud de eslabon y lc=distancia al centro de masa del eslabón
syms pi g a cero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Velocidades lineales y angulares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creamos el vector de coordenadas articulares
Q= [th1];
%Creamos el vector de velocidades articulares
Qp= [th1p];
%Creamos el vector de aceleraciones articulares
Qpp= [th1pp];
%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP=[0];

%Número de grado de libertad del robot
GDL= size(RP,2);
GDL_str= num2str(GDL);

%Articulación 1 
%Posición de la articulación 1 respecto a 0
P(:,:,1)= [l1*cos(th1);l1*sin(th1);0];
%Matriz de rotación de la junta 1 respecto a 0
R(:,:,1)= [cos(th1) -sin(th1)  0;
           sin(th1)  cos(th1)  0;
           0         0         1];

%Creamos un vector de ceros
Vector_Zeros= zeros(1, 3);

%Inicializamos las matrices de transformación Homogénea locales
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las matrices de transformación Homogénea globales
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las posiciones vistas desde el marco de referencia inercial
PO(:,:,GDL)= P(:,:,GDL); 
%Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
RO(:,:,GDL)= R(:,:,GDL); 


%Calculamos las matrices de transformación homogeneas globales y locales
for i = 1:GDL
    i_str= num2str(i);
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);

   %Globales
    try
       T(:,:,i)= T(:,:,i-1)*A(:,:,i);
    catch
       T(:,:,i)= A(:,:,i);
    end
    T(:,:,i)= simplify(T(:,:,i));
    RO(:,:,i)= T(1:3,1:3,i);
    PO(:,:,i)= T(1:3,4,i);
end

%Calculamos el jacobiano lineal de forma analítica
Jv_a(:,GDL)=PO(:,:,GDL);
Jw_a(:,GDL)=PO(:,:,GDL);

for k= 1:GDL
    if RP(k)==0 
       %Para las juntas de revolución
        try
            Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
            Jw_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)= cross([0,0,1], PO(:,:,GDL));%Matriz de rotación de 0 con respecto a 0 es la Matriz Identidad, la posición previa tambien será 0
            Jw_a(:,k)=[0,0,1];%Si no hay matriz de rotación previa se obtiene la Matriz identidad
         end
    else
        %Para las juntas prismáticas
        try
            Jv_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)=[0,0,1];
        end
            Jw_a(:,k)=[0,0,0];
     end
 end    

%Obtenemos SubMatrices de Jacobianos
Jv_a= simplify (Jv_a);
Jw_a= simplify (Jw_a);

%Matriz de Jacobiano Completa
%disp('Matriz de Jacobiano');
Jac= [Jv_a;
      Jw_a];
Jacobiano= simplify(Jac);

%Obtenemos vectores de Velocidades Lineales y Angulares

%Velocidad lineal
disp('Velocidad lineal obtenida mediante el Jacobiano lineal');
V=simplify (Jv_a*Qp);
pretty(V);

%Velocidad angular
disp('Velocidad angular obtenida mediante el Jacobiano angular');
W=simplify (Jw_a*Qp);
pretty(W);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Energía Cinética
%%%%%%%%%%%%%%%%%%%%%%%%%%Omitimos la división de cada lc%%%%%%%%%%%%%%%

%Distancia del origen del eslabón a su centro de masa
%Vectores de posición respecto al centro de masa
 P01=subs(P(:,:,1)/2, l1, lc1);%La función subs sustituye l1 por lc1 en 
                             %la expresión P(:,:,1)/2
                             
%Creamos matrices de inercia para cada eslabón
I1=[Ixx1 0 0; 
    0 Iyy1 0; 
    0 0 Izz1];

%Función de energía cinética
%Extraemos las velocidades lineales en cada eje
V=V(t); Vx= V(1,1); Vy= V(2,1); Vz= V(3,1);

%Extraemos las velocidades angular en cada ángulo de Euler
W=W(t); W_pitch= W(1,1); W_roll= W(2,1); W_yaw= W(3,1);

%Calculamos las velocidades para el eslabón1 %

%%No es necesario realizar el calculo nuevametne ya que nuestro 
%%robot cuenta con un solo eslabón. 

%Calculamos la energía cinética para cada el eslabón %%%%%%

%Eslabón 1
V1_Total= V+cross(W,P01);
K1= (1/2*m1*(V1_Total))'*((V1_Total)) + (1/2*W)'*(I1*W);
%disp('Energía Cinética en el Eslabón 3');
K1= simplify (K1);
%pretty (K3);

disp('Energía cinética total: ');
K_Total= simplify (K1);
pretty (K_Total);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Energía Potencial total y Lagrangiano 
%%%%%%%%%%%%%%%%%%%%%%%%%%Omitimos la división de cada lc%%%%%%%%%%%%%%%

%Obtenemos las alturas respecto a la gravedad
h1= P01(2); %Tomo la altura paralela al eje y

U1=m1*g*h1;

%Calculamos la energía potencial total
U_Total= U1;

disp('Energia potencial total: '); 
pretty(U_Total); 

disp('Modelo de energia: ')
H = simplify(K_Total+U_Total); 
pretty(H)

%Obtenemos el Lagrangiano
Lagrangiano= simplify (K_Total-U_Total);
disp('Lagrangiano: ')
pretty (Lagrangiano);

%%%%%%%%%%%%%%%%%%%%%%Ecuaciones de Movimiento%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lagrangiano derivado con respecto a la primera coordenada generalizada de
%velocidad

%Definimos un vector columna de derivadas con respecto al tiempo
%En este vector agrego las velocidades y aceleraciones
%Derivadas respecto al tiempo

Qd=[th1p(t); th1pp(t)];

%Obtenemos las derivadas de la velocidad en la primera coordenada
%generalizada

dQ1=[diff(diff(Lagrangiano,th1p), th1),... %Derivamos con respecto a la primera velocidad generalizada th1p para las 3 posiciones articulaciones
    diff(diff(Lagrangiano,th1p), th1p)];%Derivamos con respecto a la primera velocidad generalizada th1p para las 3 velocidades articulaciones

%Definimos el torque 1
t1= dQ1*Qd- diff(Lagrangiano, th1);

%Generación del Modelo Dinámico en forma matricial

%Matriz de Inercia

%Extraemos coeficientes de aceleraciones

M=[diff(t1, th1pp)];
rank(M);
M=M(t);

%Fuerzas Centrípetas y de Coriolis
 
%Definimos Mp
Mp= diff(M, th1p)*Qp; 

%Definimos la energía cinética en su forma matricial
k=1/2*transpose(Qp)*M*Qp;

%Definimos dk
dk=[diff(k, th1)];

%Fuerzas centrípetas y de Coriolis
C= Mp*Qp-dk;

%Par Gravitacional
%se sustituyen las velocidades y aceleraciones por 0
r=cero;
a1=subs(t1, th1p, r);

%Torque gravitacional en el motor 1
G1=a1;

%Par gravitacional
G=G1; 

disp('Matriz de inercia: ')
pretty(M); 

disp('Fuerzas centripetas y coriolis: ')
pretty(C); 

disp('Par gravitacional')
pretty(G); 
