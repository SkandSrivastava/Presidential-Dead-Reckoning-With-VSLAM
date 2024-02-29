function [Pos_G,Vel_G,Att_G]=INS(Acc,Gyr)

%=========================================================================
% Calcula algoritmo INS integrando datos de IMU 
%    Asume: - Posicion y orientaci�n inicial conocida
%           
%    Ejemplo uso: [Pos_G,Vel_G,Att_G]=INS(Acc,Gyr);
%
% INPUT:
%     Acc:       Matriz [nx3] de "n" muestras de las aceleraciones medidas con los acelerometros en X,Y y Z (m/s^2)
%     Gyr:       Matriz [nx3] de "n" muestras de las velocidades angulares medidas con giroscopos en X,Y y Z  (rad/s)
% OUTPUT:
%     Pos_G:     Posicion verdadera de IMU en sistema navegacion global G;   Po [nx3] X,Y,Z [m]
%     Vel_G:     Velocidad verdadera de IMU en sistema navegacion global G;  Ve [nx3] Vx,Vy,Vz [m/s]
%     Att_G:     Orientaci�n verdadera de IMU en frame global G;            Att [nx3] Attitude en roll,pitch y yaw [�]



%===============================================================================================================
% reservar espacio para variables
num_muestras=size(Acc,1);
tiempo_experimento=Acc(end,4)-Acc(1,4);
fs=floor(num_muestras/tiempo_experimento);  % samples/s o Hz

Rot_GS=zeros(3,3,num_muestras);
A=zeros(3,3,num_muestras);
Acc_G=zeros(num_muestras,3);
Acc_G_sin_g=zeros(num_muestras,3);
Vel_G=zeros(num_muestras,3);
Pos_G=zeros(num_muestras,3);

% Inicializar 
w=1:ceil(20*fs);  % ventana asumo al menos 20 segundos de reposo absoluto
Pos_G_ini=[0 0 0];
gravedad=mean(sqrt(Acc(w,1).^2+Acc(w,2).^2+Acc(w,3).^2));  % m/s^2
% Iniciliazar Matriz Rotaci�n:
Acc_mean=[mean(Acc(w,1)),mean(Acc(w,2)),mean(Acc(w,3))];
Att_G_ini(1)=atan2(Acc_mean(2),Acc_mean(3))*180/pi;  % Roll en grados
Att_G_ini(2)=-atan2(Acc_mean(1),sqrt(Acc_mean(2)^2+Acc_mean(3)^2))*180/pi;  % Pitch en grados
Att_G_ini(3)=0;

% matriz de rotacion (DCM) inicial
Rot_z=[cos(Att_G_ini(3)*pi/180) -sin(Att_G_ini(3)*pi/180) 0;  sin(Att_G_ini(3)*pi/180) cos(Att_G_ini(3)*pi/180) 0; 0 0 1];% from S-Sensor to G-Global
Rot_y=[cos(Att_G_ini(2)*pi/180) 0 sin(Att_G_ini(2)*pi/180); 0 1 0; -sin(Att_G_ini(2)*pi/180)   0 cos(Att_G_ini(2)*pi/180)];% from S-Sensor to G-Global
Rot_x=[1 0 0; 0 cos(Att_G_ini(1)*pi/180) -sin(Att_G_ini(1)*pi/180); 0  sin(Att_G_ini(1)*pi/180) cos(Att_G_ini(1)*pi/180)];% from S-Sensor to G-Global
Rot_GS(:,:,1)=Rot_z(:,:)*Rot_y(:,:)*Rot_x(:,:);  % from S-Sensor to G-Global BUENA

%--------------------------------------------------------------------------
%  INS
aprox=6;  % Uso este m�todo de integraci�n del giro (el mas exacto)
for i=1:num_muestras
    % Actualizar matriz de orientaci�n con valores de los "Gyros"
    SkewGyros=[0 -Gyr(i,3)  Gyr(i,2);...
        Gyr(i,3)  0  -Gyr(i,1);...
        -Gyr(i,2) Gyr(i,1) 0]; % rad/s
    if i>1
        switch aprox
            case 1
                Rot_GS(:,:,i)=Rot_GS(:,:,i-1)*(eye(3) + SkewGyros/fs ); %Pade(1,0)
            case 2
                Rot_GS(:,:,i)=Rot_GS(:,:,i-1)*(2*eye(3)+SkewGyros/fs)*inv(2*eye(3)-SkewGyros/fs); % Pade(1,1)
            case 3
                Rot_GS(:,:,i)=Rot_GS(:,:,i-1)*(3*eye(3)+2*SkewGyros/fs+0.5*(SkewGyros/fs)^2)*inv(3*eye(3)-SkewGyros/fs); % Pade(2,1)
            case 4
                Rot_GS(:,:,i)=Rot_GS(:,:,i-1)*(2*eye(3)+SkewGyros/fs+(1/6)*(SkewGyros/fs)^2)*inv(2*eye(3)-SkewGyros/fs+(1/6)*(SkewGyros/fs)^2); % Pade(2,2)
            case 5
                Rot_GS(:,:,i)=Rot_GS(:,:,i-1)*(2*eye(3)+SkewGyros/fs+(1/5)*(SkewGyros/fs)^2+(1/60)*(SkewGyros/fs)^3)*inv(2*eye(3)-SkewGyros/fs+(1/5)*(SkewGyros/fs)^2-(1/60)*(SkewGyros/fs)^3); % Pade(3,3)
            case 6
                Rot_GS(:,:,i)=Rot_GS(:,:,i-1)*expm(SkewGyros/fs); %Funcion "expm" que es un Pade(6,6) BUENA
            case 7
                sig=(sqrt(Gyr(i,:)*Gyr(i,:)')/fs);
                if sig==0
                    Rot_GS(:,:,i)=Rot_GS(:,:,i-1);
                else
                    Rot_GS(:,:,i)=Rot_GS(:,:,i-1)*(eye(3)+sin(sig)/(sig)*SkewGyros/fs+(1-cos(sig))/(sig)^2*SkewGyros*SkewGyros/fs/fs); % sin cos aprox
                end
        end
    end
    % Roll, Ptich and Yaw
    Att_G(i,1)=(atan2(Rot_GS(3,2,i),Rot_GS(3,3,i)))*180/pi;  % Roll en grados
    Att_G(i,2)=(-asin(Rot_GS(3,1,i)))*180/pi;  % Pitch
    Att_G(i,3)=(atan2(Rot_GS(2,1,i),Rot_GS(1,1,i)))*180/pi;  %Yaw
    
    %.........................................................
    % Transformar medidas sensoriales de S a G
    Acc_G(i,1:3)=(Rot_GS(:,:,i)*Acc(i,1:3)')';  %(en G: Noth-West-Up)
    
    %.........................................................
    % Eliminar gravedad terrestre (que se lee en aceler�metros)
    Acc_G_sin_g(i,:)=Acc_G(i,:) - [0 0 gravedad];  %(en G: Noth-West-Up)
    % se suma 9.8 m/s^2 si G es {north, east, down}, o se resta si G es {north, west, up}
    
    %.........................................................
    % Integrar doblemente aceleraciones (ahora ya sin la gravedad)
    if i==1 
        Vel_G(1,:)=[0 0 0];
        Pos_G(1,:)=Pos_G_ini;
    else
        % Velocidad:
        %    Vel_G(i,:)=Vel_G(i-1,:) + Acc_G_sin_g(i,:)/fs; % integracion simple
        Vel_G(i,:)=Vel_G(i-1,:) + 0.5*(Acc_G_sin_g(i,:)+Acc_G_sin_g(i-1,:))/fs;  % integracion trapezoidal
        %   Vel_G(i,:)=Vel_G(i-1,:) + (1/6)*(Acc_G_sin_g(i,:)+4*Acc_G_sin_g(i-1,:)+Acc_G_sin_g(i-2,:))/fs; % integracion Simpsons's rule
                
        %Posicion:
        %  Pos_G(i,:)=Pos_G(i-1,:) + Vel_G(i,:)/fs; % integracion simple
        Pos_G(i,:)=Pos_G(i-1,:) + 0.5*(Vel_G(i,:)+Vel_G(i-1,:))/fs;  % integracion trapezoidal
        %  Pos_G(i,:)=Pos_G(i-1,:) + (1/6)*(Vel_G(i,:)+4*Vel_G(i-1,:)+Vel_G(i-2,:))/fs; % integracion Simpsons's rule
    end
end

% TEXTO para explicar en clase sobre el INS:
% %--------------------------------------------------------------------------
% %  INS
% %  Acc: IMU Acelerometer signal (3 axis) [m/s^2]
% %  Gyr: IMU Angular Rate signal (3 axis) [rad/s]
% %  fs:  Frecuencia muestreo [Hz]
%
% % Actualizar matriz de orientaci�n con valores de los "Gyros"
% SkewGyros=[0 -Gyr(i,3)  Gyr(i,2);...
%         Gyr(i,3)  0  -Gyr(i,1);...
%         -Gyr(i,2) Gyr(i,1) 0]*pi/180; % rad/s
% Rot_nb(:,:,i)=Rot_nb(:,:,i-1)*(2*eye(3)+SkewGyros/fs)*inv(2*eye(3)-SkewGyros/fs); % Pade(1,1)
%
% %.........................................................
% % Transformar medidas sensoriales de "b-frame" a "n-frame"
% Acc_n(i,:)=(Rot_nb(:,:,i)*Acc(i,:)')';  %(en "n-frame": Noth-West-Up)
%
% %.........................................................
% % Eliminar gravedad terrestre (que se lee en aceler�metros)
% Acc_n_sin_g(i,:)=Acc_n(i,:) - [0 0 gravedad];  %(en "n-frame": Noth-West-Up)
%
% %.........................................................
% % Integrar doblemente aceleraciones (ahora ya sin la gravedad)
% Vel_n(i,:)=Vel_n(i-1,:) + 0.5*(Acc_n_sin_g(i,:)+Acc_n_sin_g(i-1,:))/fs;  % integracion trapezoidal
% Pos_n(i,:)=Pos_n(i-1,:) + 0.5*(Vel_n(i,:)+Vel_n(i-1,:))/fs;  % integracion trapezoidal



