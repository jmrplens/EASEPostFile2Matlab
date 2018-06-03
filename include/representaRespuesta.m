function representaRespuesta(handles)

axes(handles.plotfig)
Resp = handles.SPLm{handles.Index,13};
Frecuencia = [125,250,500,1000,2000,4000,8000];
Fuente = string(handles.SPLm(handles.Index,2));
% Pasar de tercios de octava a octavas
Octavas = zeros(1,7);
a = 0;
for i=2:3:21 
    a = a+1;
    Octavas(a) = 10*log10(sum(10.^(Resp(i-1)/10)+10.^(Resp(i)/10)+10.^(Resp(i+1)/10)));  
end

bar(Octavas)
set(gca,'XTickLabel', Frecuencia);

title(sprintf('Respuesta en frecuencia de la fuente %s',Fuente))
xlabel('Frecuencia (Hz)')
ylabel('Nivel SPL (dB)')
end