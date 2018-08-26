function importarPostFiles(hObject,handles)
% Almacena en la aplicacion (handles.SPLm) la variables
% Es una variable tipo cell que en cada fila almacena los datos extraidos
% de cada fila. Si hay 20 archivos, el cell 'SPLm' tendr· 20 filas.
% En las ultimas lineas se puede observar que contiene cada columna.
%
% Distancia -> Vector de distancias entre cada receptor y la fuente
% PosReceptor -> Vector con la posicion del receptor (X,Y,Z)
% PosFuente -> Vector con la posicion de la fuente (X,Y,Z)
% Fuente -> Vector con la ID de fuente recibida cada receptor
% Receptor -> Vector con las ID de cada receptor
% Volumen -> Variable que contiene el valor del volumen de la sala
% Superficie -> Variable que contiene el valor de la superficie
% Temperatura -> Variable que contiene el valor de la temperatura
% Humedad -> Variable que contiene el valor de la humedad
% SPLm -> Cell que almacena la historia temporal de cada receptor
% RayoDirecto -> Tiempo de llegada del primer rayo a cada receptor
% RetFuente ->  Retardo de la fuente
% NumRayos -> Numero de rayos recibidos por cada receptor
% Reverberacion -> Tiempo de reverberacion de cada archivo
% Absorcion -> Coeficiente de absorcion de cada archivo
% Directividad -> Directividad de las fuentes
% RuidoFondo -> Ruido de fondo de cada archivo
% SPLFuente -> Nivel SPL de las fuentes

%% Importacion de datos
ListaArchivos = handles.ListaArchivos;

%multiWaitbar(handles,'Importando archivos',N/size(PARAMtxt,1),'color','b');
% Inicializacion de variables
Distancia = zeros(size(ListaArchivos,2),1);
PosReceptor = zeros(size(ListaArchivos,2),3);
PosFuente = zeros(size(ListaArchivos,2),3);
RayoDirecto = zeros(size(ListaArchivos,2),1);
RetFuente = zeros(size(ListaArchivos,2),1);
NumRayos = zeros(size(ListaArchivos,2),1);
Reverberacion = zeros(size(ListaArchivos,2),21);
Absorcion = zeros(size(ListaArchivos,2),21);
Directividad = zeros(size(ListaArchivos,2),21);
RuidoFondo = zeros(size(ListaArchivos,2),21);
SPLFuente = zeros(size(ListaArchivos,2),21);
handles.SPLm = cell(size(ListaArchivos,2),21);
for i=1:size(ListaArchivos,2)
    multiWaitbar(handles,'Importando datos',i/size(ListaArchivos,2),'color','b');
    % Distancia entre receptores y fuente
    A = EASEPost2Distan(ListaArchivos{i});
    Distancia(i) = sqrt(...
        (A(1,1)-A(2,1))^2 + (A(1,2)-A(2,2))^2 + (A(1,3)-A(2,3))^2);
    PosReceptor(i,:) = A(1,:);
    PosFuente(i,:) = A(2,:);
    
    % Nombre de la fuente y de los receptores
    A = EASEPost2Info(ListaArchivos{i});
    Fuente{i} = char(A(2));
    Receptor{i} = char(A(1));
    
    % Parametros de la sala
    % ToDo: Almacenar la primera vez y el resto comprobar que son iguales, si
    % no son iguales mostrar error y detener la importaciÛn
    A = EASEPost2Sala(ListaArchivos{i});
    Volumen = A(1);
    Superficie = A(2);
    Temperatura = A(3);
    Humedad = A(4);
    
    % Datos varios
    A = EASEPost2Otros(ListaArchivos{i});
    RayoDirecto(i) = A(1);
    RetFuente(i) = A(2);
    NumRayos(i) = A(3);
    
    % RT, Absorcion de la sala, Directividad fuente
    A = EASEPost2Param(ListaArchivos{i});
    Reverberacion(i,:) = A(1,:);
    Absorcion(i,:) = A(2,:);
    Directividad(i,:) = A(3,:);
    
    % Ruido de fondo, SPL fuente
    A = EASEPost2Param2(ListaArchivos{i});
    RuidoFondo(i,:) = A(1,:);
    SPLFuente(i,:) = A(2,:);
    
    % Obtener historia temporal para cada receptor
    % Obtiene una matriz con la columna de delay y los niveles en tercios de
    % octava
    A = EASEPost2Matrix(ListaArchivos{i});
    % Eliminar el delay de llegada, restando el valor minimo de llegada a todo
    % el vector de tiempos. Asi aseguramos que la integracion por rangos
    % temporales sea correcta.
    Tnew = A(:,1) - RayoDirecto(i);
    % Almacenar en el cell 'SPLm' la informacion del archivo actual 'i'
    handles.SPLm{i,1} = []; % Aqui ira el SPL por milisegundo
    handles.SPLm{i,2} = Fuente{i};
    handles.SPLm{i,3} = Receptor{i};
    handles.SPLm{i,4} = PosFuente(i,:);
    handles.SPLm{i,5} = PosReceptor(i,:);
    handles.SPLm{i,6} = Distancia(i);
    handles.SPLm{i,7} = RetFuente(i);
    handles.SPLm{i,8} = NumRayos(i);
    handles.SPLm{i,9} = RayoDirecto(i);
    handles.SPLm{i,10} = Reverberacion(i,:);
    handles.SPLm{i,11} = Absorcion(i,:);
    handles.SPLm{i,12} = RuidoFondo(i,:);
    handles.SPLm{i,13} = SPLFuente(i,:);
    handles.SPLm{i,14} = Directividad(i,:);
    handles.SPLm{i,15} = Volumen;
    handles.SPLm{i,16} = Superficie;
    handles.SPLm{i,17} = Temperatura;
    handles.SPLm{i,18} = Humedad;
    % Historia temporal completa
    handles.SPLm{i,19} = [Tnew,A(:,2:end)];
    % Historia temporal procesada por milisegundo en bandas de octava
    handles.SPLm{i,20} = [];
    % Historia temporal procesada por milisegundo en tercios de octava
    handles.SPLm{i,21} = [];
    
end

multiWaitbar(handles,'Importando datos','Close');
% En primer lugar agrupa los niveles en bloques de 1 milisegundo, ya que hay
% varios niveles por milisegundo, y promedia los valores de cada bloque
for i=1:numel(handles.SPLm(:,1))
    
    multiWaitbar(handles,'Procesando historias temporales',i/numel(handles.SPLm(:,1)),'color','r');
    for j=1:max(ceil(handles.SPLm{i,19}(:,1)))
        Ini = find(handles.SPLm{i,19}(:,1)>=j-1,1); % Inicio del bloque de 'j' milisegundos
        Fin = find(handles.SPLm{i,19}(:,1)>j,1)-1; % Fin del bloque de 'j' milisegundos
        % Si ha llegado al final del vector, fin estara vacio y inicio y fin
        % seran iguales para obtener el valor ultimo del vector
        if isempty(find(handles.SPLm{i,19}(:,1)>j,1)); Fin=Ini; end
        % Si el bloque del milisegundo 'j' existe, ini sera menor o igual a fin,
        % si no existe ini sera mayor a fin y no habra valor en ese milisegundo.
        if Ini<=Fin
            Valores = handles.SPLm{i,19}(Ini:Fin,2:end);
            
            if max(Valores(:)) > 0
            % Suma el bloque por octavas y vuelve a sumar para obtener un valor
            ValMS(j) = 10*log10(sum(sum(10.^(Valores/10),1)));
            % Suma los rayos recibidos en el milisegundo 'j', y lo pasa de
            % tercios de octava a bandas de octava
            ValMSthirdoct(j,:) = sum(10.^(Valores/10),1);
            ValMSoct(j,:) = 10*log10([sum(ValMSthirdoct(j,1:3)),... % 125 Hz
                sum(ValMSthirdoct(j,4:6)),...  % 250 Hz
                sum(ValMSthirdoct(j,7:9)),...  % 500 Hz
                sum(ValMSthirdoct(j,10:12)),...% 1 kHz
                sum(ValMSthirdoct(j,13:15)),...% 2 kHz
                sum(ValMSthirdoct(j,16:18)),...% 4 kHz
                sum(ValMSthirdoct(j,19:21))]); % 8 kHz
            ValMSthirdoct(j,:) = 10*log10(ValMSthirdoct(j,:));
            else
              ValMS(j) = 0;
              ValMSoct(j,:) = 0;
              ValMSthirdoct(j,:) = 0;
            end
        else
            ValMS(j)=NaN;
            ValMSoct(j,:) = NaN;
            ValMSthirdoct(j,:) = NaN;
        end
        
    end
    % Todos los valores 0 y menores a 0 se dejan vacios
    %ValMS(ValMS==0) = NaN;
    %ValMSoct(ValMSoct==0) = NaN;
    % Si se utiliza matlab 2014b o anterior realiza la interpolacion
    % manualmente, sino utiliza la funcion fillmising
    if verLessThan('matlab','8.7')
        % Interpolacion manual global
        bd=isnan(ValMS);
        gd=find(~bd);
        bd([1:(min(gd)-1) (max(gd)+1):end])=0;
        ValMS(bd)=interp1(gd,ValMS(gd),find(bd));
        % Interpolacion manual por bandas
        for jj=1:7
            bd=isnan(ValMSoct(:,jj));
            gd=find(~bd);
            bd([1:(min(gd)-1) (max(gd)+1):end])=0;
            ValMSoct(bd,jj)=interp1(gd,ValMSoct(gd,jj),find(bd));
        end
        for jj=1:21
            bd=isnan(ValMSthirdoct(:,jj));
            gd=find(~bd);
            bd([1:(min(gd)-1) (max(gd)+1):end])=0;
            ValMSthirdoct(bd,jj)=interp1(gd,ValMSthirdoct(gd,jj),find(bd));
        end
    else
        % Interpolar valores en los milisegundos que haya valor 0, desde el inicio
        % al final del vector
        ValMS = fillmissing(ValMS,'linear','EndValues','extrap');
        ValMSoct = fillmissing(ValMSoct,'linear','EndValues','extrap');
        ValMSthirdoct = fillmissing(ValMSthirdoct,'linear','EndValues','extrap');
    end
    handles.SPLm{i,1} = ValMS';
    handles.SPLm{i,20} = ValMSoct;
    handles.SPLm{i,21} = ValMSthirdoct;
end
multiWaitbar(handles,'Procesando historias temporales','Close');
% Obtener el tiempo maximo valor elegible para representacion. Obtiene el
% maximo de cada receptor y despues busca el minimo entre estos.
handles.MaxElegible = min(cell2mat(cellfun( @(val) numel(val), handles.SPLm(:,1),'UniformOutput',0)))-1;

% Ordenar los datos segun la ID de fuente y despues por receptor. Asi se
% tiene agrupados los datos segun la fuente, necesario para elegir una
% fuente u otra en el programa.
handles.SPLm = sortrows(handles.SPLm,[2,3]);
% Actualizar variable en la aplicacion
guidata(hObject, handles);

end

