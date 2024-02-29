function [Num_steps,StanceBegins_idx,StancePhase,idx_fig]=StepDetection_Acel_smartphone(Acc,ver,idx_fig)

%=========================================================================
% Algorithm que usa solo los datos de los acelerómetros. Basado en: [Stirling2003]
% Pasos: 1)Acc magnitude + 2) LowPass + 3)Local Variance + 4) "Begin-of-Stance" event: Detection of sudden variance drop.
%
% INPUT:
%      Acc:      Matriz [nx4], con "n" muestras de los Acelerometros para los ejes X, Y y Z, time en las columnas 1, 2, 3 y 4, respectivamente.
%      ver:      0: no ver figuras;  distinto de cero: verlas
%      idx_fig:  entero que indica en que figura hay que hacer la primera representación
% OUTPUT:
%     Num_steps: entero indicando el número de pasos detectados
%     StanceBegins_idx:  Events when a step occurs
%     StancePhase:       1 or 0. It is 1 if on stance, 0 otherwise
%     idx_fig:   entero que indica en que figura se debería pintar la proxima vez


%..............preprocesado..................
num_samples=size(Acc,1);
tiempo_experimento=Acc(end,4)-Acc(1,4);
freq_Acc=ceil(num_samples/tiempo_experimento);  % samples/s o Hz
Acc_mag=sqrt(Acc(:,1).^2+Acc(:,2).^2+Acc(:,3).^2);
orden_filtro=4;
[b,a] = butter(orden_filtro,2/(freq_Acc/2),'low'); % diseño filtro IIR paso baja de orden 4 (freq corte 2 Hz)
[Acc_mag_filt,zf]=filter(b,a,Acc_mag,[[9.73783786570634;-22.8430473945528;18.4412646506037;-5.04055548881593]]);

% .............Detectar pasos ..................
umbral_Acc=0.4; % umbral de 0.4 m/s^2
umbral_Acc_descarte=2.0; % umbral superior a 1.5 m/s^2 que indica meneo excesivo
gravity=9.8;
Acc_filt_binary=zeros(1,length(Acc));
Acc_filt_detrend=zeros(1,length(Acc));
for ii=2:length(Acc)
    gravity=0.999*gravity+0.001*Acc_mag(ii);  % calculo gravedad experimentalmente (pues en S3 me da 9.7 en vez de lo esperado que es 9.8)
    Acc_filt_detrend(ii)=Acc_mag_filt(ii)-gravity;
    if Acc_filt_detrend(ii)>umbral_Acc && Acc_filt_detrend(ii)<umbral_Acc_descarte
        Acc_filt_binary(ii)=1; % fases de subida del cuerpo (inicio paso), pongo "1"
    else
        if Acc_filt_detrend(ii)<-umbral_Acc
            if Acc_filt_binary(ii-1)==1 % si despues de un "1" no hubiese un valor intermedio para generar el "0" entonces se lo pongo de manera forzada
                Acc_filt_binary(ii)=0;
            else  % caso normal: ya he puesto el "1", luego en "0", y ahora el "-1"
                Acc_filt_binary(ii)=-1; % fases de bajada del cuerpo (fin paso)
            end
        else
            Acc_filt_binary(ii)=0; % si está entre el umbral superior y el inferior => lo marco con "0"
        end
    end
    
end
stepcount=0;
StanceBegins_idx=[];
StepDect=zeros(1,length(Acc));
steps=zeros(1,length(Acc))*NaN;
window=ceil(0.4*freq_Acc); % muestras en ventana para considerar 0.4 segundos
for ii=(window+2):length(Acc)
    % detecto pasos si:
    if ( Acc_filt_binary(ii)==-1  && Acc_filt_binary(ii-1)==0  && (sum(Acc_filt_binary(ii-window:ii-2))>1 ))
        StepDect(ii)=1;
        stepcount=stepcount+1; % contador de pasos
        StanceBegins_idx(stepcount)=ii;
    end
    if StepDect(ii)
        steps(ii)=0;
    else
        steps(ii)=NaN;
    end
end

% Todas las muetras de apoyo
StancePhase=zeros(num_samples,1);
for ii=StanceBegins_idx
   StancePhase(ii:ii+9,1)=ones(10,1);  % asumo que fase de apoyo son las 10 siguientes muestras al comienzo de fase de apoyo (StanceBegins_idx) 
end

% Cálculo parámetros:
Num_steps=length(StanceBegins_idx); % numero de pasos contados

figure(idx_fig); idx_fig=idx_fig+1;
vect_smaples=1:num_samples;
plot(vect_smaples,Acc_mag,'r-'); hold on;  title('"SL+theta PDR": Acclerometer processing for Step detection'); ylabel('Accleration'); xlabel('time (seconds)');
plot(vect_smaples,Acc_mag_filt,'b-');
plot(vect_smaples,Acc_filt_detrend,'c-');
plot(vect_smaples,Acc_filt_binary,'gx-');
plot(vect_smaples,steps,'ro','MarkerSize',8); hold off;
legend({'|Acc|','lowpass(|Acc|)','detrend(lowpass(|Acc|))','Binary','detected Steps'});

