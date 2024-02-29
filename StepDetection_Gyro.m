function [Num_steps,StanceSample,StancePhase,idx_fig]=StepDetection_Gyro(Gyr,ver,idx_fig)

%=========================================================================
%  Algorithm 1: Step Detetion using Gyros
%  Rate gyro magnitude + Umbral + Mediana + "Begin-of-Stance" event
%  Basado en: [Feliz 2009]
%
% INPUT:
%      Gyr:      Matriz [nx3], con "n" muestras de los giroscopos, para los ejes X, Y y Z, en las filas 1, 2 y 3, respectivamente.
%      ver:      0: no ver figuras;  distinto de cero: verlas
%      idx_fig:  entero que indica en que figura hay que hacer la primera representación
% OUTPUT:
%     Num_steps: entero indicando el número de pasos detectados
%     idx_fig:   entero que indica en que figura se debería pintar la proxima vez


num_samples=size(Gyr,1);

% inicializar variables
M_Gyr3=zeros(1,num_samples);

% Paso1: rate gyro magnitude
M_Gyr=sqrt(Gyr(:,1).^2+Gyr(:,2).^2+Gyr(:,3).^2);

% Paso 2: Thresholding
Th_M_Gyr=1; % 1 rad/s
M_Gyr2=int8(M_Gyr>Th_M_Gyr)*Th_M_Gyr;

% Paso 3: filtro mediana (ventana 0.20 seg; 20 muestras)
half_size=10;
M_Gyr3(1:half_size)=M_Gyr2(1:half_size);
for i=half_size+1:length(M_Gyr2)-half_size
    valores=M_Gyr2(i-half_size:i+half_size-1);
    M_Gyr3(i)=median(valores);
end
M_Gyr3(length(M_Gyr2)-half_size+1:length(M_Gyr2))=...
    M_Gyr2(length(M_Gyr2)-half_size+1:length(M_Gyr2));
M_Gyr4=M_Gyr3;

% % Paso 3: filtro dilatacion+erosion (ventana 0.20 seg; 20 muestras)
% half_size=6;
% % dilatacion
% M_Gyr3(1:half_size)=M_Gyr2(1:half_size);
% for i=half_size+1:length(M_Gyr2)-half_size
%     valores=M_Gyr2(i-half_size:i+half_size-1);
%     if max(valores)==Th_M_Gyr
%         M_Gyr3(i)=Th_M_Gyr;
%     else
%         M_Gyr3(i)=M_Gyr2(i);
%     end
% end
% M_Gyr3(length(M_Gyr2)-half_size+1:length(M_Gyr2))=M_Gyr2(length(M_Gyr2)-half_size+1:length(M_Gyr2));
% % erosion
% M_Gyr4(1:half_size)=M_Gyr3(1:half_size);
% for i=half_size+1:length(M_Gyr2)-half_size
%     valores=M_Gyr3(i-half_size:i+half_size-1);
%     if min(valores)==0
%         M_Gyr4(i)=0;
%         else
%         M_Gyr4(i)=M_Gyr3(i);
%     end
% end
% M_Gyr4(length(M_Gyr3)-half_size+1:length(M_Gyr3))=M_Gyr3(length(M_Gyr3)-half_size+1:length(M_Gyr3));

% Paso 4: Begin of stance event
M_Gyr5(1)=0;
M_Gyr5(2:length(M_Gyr4))=M_Gyr4(2:length(M_Gyr4))<M_Gyr4(1:length(M_Gyr4)-1);
StanceSample=find(M_Gyr5);

% Todas las muetras de apoyo
StancePhase=zeros(num_samples,1);
for i=StanceSample
   StancePhase(i:i+9,1)=ones(10,1);  % asumo que fase de apoyo son las 10 siguientes muestras al comienzo de fase de apoyo (StanceBegins_idx) 
end


% Calculo parámetros:
Num_steps=length(StanceSample); % numero de pasos contados

if ver
    % ver magnitudes del giroscopo
    figure(idx_fig);
    plot(Gyr(:,1),'r'); hold on;
    plot(Gyr(:,2),'g');
    plot(Gyr(:,3),'b');
    plot(M_Gyr,'k','LineWidth',2);
    ylabel('rad/s'); title('Angular rate');
    xlabel('samples');
    legend('w_x','w_y','w_z','Magnitude');
    idx_fig=idx_fig+1;hold off;
    
    % ver etapas intermedias de deteccion de pasos
    figure(idx_fig);
    plot(M_Gyr,'k','LineWidth',1); hold on;
    plot(M_Gyr2,'g','LineWidth',2);
    plot(M_Gyr4,'b','LineWidth',2);
    plot(StanceSample,zeros(1,Num_steps),'ro','MarkerSize',10,'MarkerFaceColor',[1 0 0]);
    ylabel('rad/s'); title('Step detection with angular rate thresholding');
    xlabel('samples');
    legend('Magnitude','Threshold','Median filter','Stance-Begins');
    idx_fig=idx_fig+1; hold off;
end