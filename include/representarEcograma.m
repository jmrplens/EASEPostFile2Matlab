function representarEcograma(handles)

axes(handles.plotfig)
Index = handles.Index;
SPLm = handles.SPLm(:,1);
SPLmTot = handles.SPLm(:,19);
Distancia = cell2mat(handles.SPLm(:,6))';
Fuente = string(handles.SPLm(:,2)');
Receptor = string(handles.SPLm(:,3)');


% Obtener ecograma del receptor Index, sumando las octavas
Ecogra = 10*log10(sum(10.^((SPLmTot{Index}(:,2:end)')/10)));
tiempo = SPLmTot{Index}(:,1);



%% Representacion grafica


% Ecograma
stem(tiempo,Ecogra,'Marker','none')

title({sprintf('Ecograma del receptor: %s',Receptor(Index));...
    ['\color{blue} ','Fuente: ',sprintf(': %s',Fuente(1))]})
xlabel('Tiempo (ms)')
ylabel('Nivel SPL (dB)')
% Limite del eje x con margen a ambos lados
xlim([-max(tiempo)/10,max(tiempo)+max(tiempo)/10])
end