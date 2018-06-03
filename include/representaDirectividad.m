function representaDirectividad(handles)

axes(handles.plotfig)
Resp = handles.SPLm{handles.Index,14};
Frecuencia = [125,250,500,1000,2000,4000,8000];
Fuente = string(handles.SPLm(handles.Index,2));
% Obtener solo el valor de las octavas
Octavas = zeros(1,7);
a = 0;
for i=2:3:21
    a = a+1;
    Octavas(a) = Resp(i);  
end

bar(Octavas)
set(gca,'XTickLabel', Frecuencia);

title(sprintf('Respuesta en frecuencia de la fuente %s',Fuente))
xlabel('Frecuencia (Hz)')
ylabel('DI (dB)')
end