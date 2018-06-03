function representaAbsorcion(handles)

axes(handles.plotfig)
Absor = handles.SPLm{handles.Index,11};
Frecuencia = [125,250,500,1000,2000,4000,8000];
Fuente = string(handles.SPLm(handles.Index,2));
% Pasar de tercios de octava a octavas
Octavas = zeros(1,7);
a = 0;
for i=2:3:21 
    a = a+1;
    Octavas(a) = (Absor(i-1)+Absor(i)+Absor(i+1))/3;  
end

bar(Octavas*100)
set(gca,'XTickLabel', Frecuencia);

title('Coeficientes de absorción medio')
xlabel('Frecuencia (Hz)')
ylabel('\alpha (%)')
end