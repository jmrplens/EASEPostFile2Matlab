function representaRT(handles)

axes(handles.plotfig)
RT = handles.SPLm{handles.Index,10};
Frecuencia = [125,250,500,1000,2000,4000,8000];
Fuente = string(handles.SPLm(handles.Index,2));
% Pasar de tercios de octava a octavas
Octavas = zeros(1,7);
a = 0;
for i=2:3:21 
    a = a+1;
    Octavas(a) = (RT(i-1)+RT(i)+RT(i+1))/3;  
end

bar(Octavas)
set(gca,'XTickLabel', Frecuencia);

title(sprintf('Tiempo de reverberación con la fuente %s',Fuente))
xlabel('Frecuencia (Hz)')
ylabel('Tiempo de reverberación (s)')
end