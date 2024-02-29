function [Num_steps,StanceBegins_idx,StancePhase,idx_fig]=StepDetection_Acel(Acc,ver,idx_fig)

%=========================================================================
% Algorithm que usa solo los datos de los acelerómetros. Basado en: [Stirling2003]
% Pasos: 1)Accel magnitude + 2) LowPass + 3)Local Variance + 4) "Begin-of-Stance" event: Detection of sudden variance drop.
%
% INPUT:
%      Acc:      Matriz [nx3], con "n" muestras de los Acelerometros para los ejes X, Y y Z, en las filas 1, 2 y 3, respectivamente.
%      ver:      0: no ver figuras;  distinto de cero: verlas
%      idx_fig:  entero que indica en que figura hay que hacer la primera representación
% OUTPUT:
%     Num_steps: entero indicando el número de pasos detectados
%     StanceBegins_idx:  Events when a step occurs
%     StancePhase:       1 or 0. It is 1 if on stance, 0 otherwise
%     idx_fig:   entero que indica en que figura se debería pintar la proxima vez

num_samples=size(Acc,1);

% inicializar variables
M_Acc_MediaLocal=zeros(1,num_samples);
M_Acc_VarLocal=zeros(1,num_samples);
M_Acc_VarLocal_Saltos=zeros(1,num_samples);

% Paso1: Accel magnitude
M_Acc=sqrt(Acc(:,1).^2+Acc(:,2).^2+Acc(:,3).^2);

% Paso 2: Local Variance
semi_ventana=25;
M_Acc_MediaLocal(1 : semi_ventana)=M_Acc(1 : semi_ventana);
M_Acc_VarLocal( 1 : semi_ventana )=zeros(1,semi_ventana);
for i=(semi_ventana+1):(length(M_Acc)-semi_ventana)
    M_Acc_MediaLocal(i)=sum(M_Acc(i-semi_ventana:i+semi_ventana))/(2*semi_ventana+1);
    M_Acc_VarLocal(i)=sum( ( M_Acc(i-semi_ventana:i+semi_ventana)-M_Acc_MediaLocal(i) ).^2 )/(2*semi_ventana+1);
end
M_Acc_MediaLocal( length(M_Acc)-semi_ventana+1 : length(M_Acc) )=M_Acc(length(M_Acc)-semi_ventana+1 : length(M_Acc));
M_Acc_VarLocal( length(M_Acc)-semi_ventana+1 : length(M_Acc) )=zeros(1,semi_ventana);
M_Acc_VarLocal=sqrt(M_Acc_VarLocal);  % en realidad voy a trabajar con la desviacion standard (para que quede mejor en los plots)

% Paso 3: Thresholding
Th_M_Acc_VarLocal1=1.6; % 2 (m/s^2)^2
Th_M_Acc_VarLocal2=0.8; % 1 (m/s^2)^2
M_Acc_VarLocal_Bool1=double(M_Acc_VarLocal>Th_M_Acc_VarLocal1)*Th_M_Acc_VarLocal1;
M_Acc_VarLocal_Bool2=double(M_Acc_VarLocal<Th_M_Acc_VarLocal2)*Th_M_Acc_VarLocal2;
M_Acc_VarLocal_Bool2(find(M_Acc_VarLocal_Bool2==0))=NaN; % para visualizar mejor

% Paso 4: Begin of stance event
ventana_futura=20; % 15 muestras son 0.15 seg (a 100 Hz)
ventana_pasada=ventana_futura;
M_Acc_VarLocal_Saltos(1)=0;
for i=1:length(M_Acc)-ventana_futura
      M_Acc_VarLocal_Saltos(i)=...
           M_Acc_VarLocal_Bool1(i+1)<M_Acc_VarLocal_Bool1(i)  &...  % bajada en umbral alto
           max(M_Acc_VarLocal_Bool2(i:i+ventana_futura))>0;         % proximo a umbral bajo
      if M_Acc_VarLocal_Saltos(i)>0 % Si salto detectado segun umbrales
          M_Acc_VarLocal_Saltos(i-ventana_pasada:i-1)=zeros(1,ventana_pasada);% eliminar cualquier otros posible salto en ventana pasada
      end
end
M_Acc_VarLocal_Saltos(length(M_Acc)-ventana_futura+1:length(M_Acc))=zeros(1,ventana_futura);
StanceBegins_idx=find(M_Acc_VarLocal_Saltos);

% Todas las muetras de apoyo
StancePhase=zeros(num_samples,1);
for i=StanceBegins_idx
   StancePhase(i:i+9,1)=ones(10,1);  % asumo que fase de apoyo son las 10 siguientes muestras al comienzo de fase de apoyo (StanceBegins_idx) 
end


% Cálculo parámetros:
Num_steps=length(StanceBegins_idx); % numero de pasos contados

if ver
    % ver magnitudes del acelerometro
    figure(idx_fig);
    plot(Acc(:,1),'r'); hold on;
    plot(Acc(:,2),'g');
    plot(Acc(:,3),'b');
    plot(M_Acc,'k','LineWidth',2);
    ylabel('m/s^2'); title('Acceleration');
    xlabel('samples');
    legend('a_x','a_y','a_z','Magnitude','Location','NorthWest');
    idx_fig=idx_fig+1;hold off;

    % ver etapas intermedias de deteccion de pasos
    figure(idx_fig);
    plot(M_Acc,'k','LineWidth',1); hold on;
    plot(M_Acc_MediaLocal,'r','LineWidth',2);
    plot(M_Acc_VarLocal,'g','LineWidth',2);
    plot(M_Acc_VarLocal_Bool1,'b','LineWidth',2);
    plot(M_Acc_VarLocal_Bool2,'m','LineWidth',2);
    plot(StanceBegins_idx,zeros(1,Num_steps),'ro','MarkerSize',10,'MarkerFaceColor',[1 0 0]);
    ylabel('m/s^2'); title('Step detection with Acceleration readings');
    xlabel('samples');
    handler_legend=legend('$a_i$','$\overline{a_i}$','$\sigma_{a_i}$','$B_1$','$B_2$','Step');
    set(handler_legend,'Interpreter','latex','Location','NorthWest');
    idx_fig=idx_fig+1; hold off;
end