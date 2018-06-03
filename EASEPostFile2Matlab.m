function varargout = EASEPostFile2Matlab(varargin)
% EASEPOSTFILE2MATLAB MATLAB code for EASEPostFile2Matlab.fig
%      EASEPOSTFILE2MATLAB, by itself, creates a new EASEPOSTFILE2MATLAB or raises the existing
%      singleton*.
%
%      H = EASEPOSTFILE2MATLAB returns the handle to a new EASEPOSTFILE2MATLAB or the handle to
%      the existing singleton*.
%
%      EASEPOSTFILE2MATLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EASEPOSTFILE2MATLAB.M with the given input arguments.
%
%      EASEPOSTFILE2MATLAB('Property','Value',...) creates a new EASEPOSTFILE2MATLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EASEPostFile2Matlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EASEPostFile2Matlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EASEPostFile2Matlab

% Last Modified by GUIDE v2.5 18-Apr-2018 01:21:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EASEPostFile2Matlab_OpeningFcn, ...
    'gui_OutputFcn',  @EASEPostFile2Matlab_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EASEPostFile2Matlab is made visible.
function EASEPostFile2Matlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EASEPostFile2Matlab (see VARARGIN)

% Choose default command line output for EASEPostFile2Matlab
handles.output = hObject;

addpath('include')

% Valor por defecto del tiempo de integracion
handles.Rango = 50; % milisegundos

% Textos de la ventana de carga
handles.LPROGRESSTITLE = 'Progreso';
handles.LTIMEDAYS = 'dias';
handles.LTIMEHOURS = 'horas';
handles.LTIMEMINS = 'mins';
handles.LTIMESECS = 'segs';
handles.LTIMESEC = 'seg';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EASEPostFile2Matlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EASEPostFile2Matlab_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in botonRT.
function botonRT_Callback(hObject, eventdata, handles)
representaRT(handles)


% --- Executes on selection change in listafuentes.
function listafuentes_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listafuentes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string',[])

% --- Executes on selection change in listareceptores.
function listareceptores_Callback(hObject, eventdata, handles)
handles.Index = get(handles.listareceptores,'Value');
% Receptor
set(handles.posicionreceptortexto,'String',...
    [num2str(handles.SPLm{handles.Index,5}(1)),'x',...
    num2str(handles.SPLm{handles.Index,5}(2)),'x',...
    num2str(handles.SPLm{handles.Index,5}(3))])
set(handles.tiempollegadatexto,'String',num2str(handles.SPLm{handles.Index,9}))
set(handles.distanciafuentetexto,'String',num2str(handles.SPLm{handles.Index,6}))
set(handles.rayosrecibidostexto,'String',num2str(handles.SPLm{handles.Index,8}))
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listareceptores_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string',[])


% --- Executes on button press in botonrespuestafuente.
function botonrespuestafuente_Callback(hObject, eventdata, handles)
representaRespuesta(handles)


% --- Executes on button press in botondirectividadfuente.
function botondirectividadfuente_Callback(hObject, eventdata, handles)
representaDirectividad(handles)


% --- Executes on button press in botonplotear.
function botonplotear_Callback(hObject, eventdata, handles)

% Si esta marcada la opcion de SPL vs Distancia
if get(handles.splVSdist,'Value')==1
    handles.Rango = str2double(get(handles.rangotemporal,'String'));
    % Muestra el mensaje de carga de representacion
    set(handles.cargandoPlot,'visible','on')
    if get(handles.doscampos,'value')==1    % Util/Perjudicial
    representarSPLvsDist2(hObject,handles)
    elseif get(handles.trescampos,'value')==1 % Directo/Early/Late
        representarSPLvsDist3(hObject,handles)
    end
    set(handles.cargandoPlot,'visible','off')
end
if get(handles.ecogramaBoton,'Value')==1
    representarEcograma(handles)
end
if get(handles.claridadboton,'Value')==1
    handles.Rango = str2double(get(handles.rangotemporal,'String'));
    set(handles.cargandoPlot,'visible','on')
    representarCt(hObject,handles)
    set(handles.cargandoPlot,'visible','off')
end
if get(handles.definicionboton,'Value')==1
    handles.Rango = str2double(get(handles.rangotemporal,'String'));
    set(handles.cargandoPlot,'visible','on')
    representarDt(hObject,handles)
    set(handles.cargandoPlot,'visible','off')
end
if get(handles.sonoridadboton,'Value')==1
    set(handles.cargandoPlot,'visible','on')
    representarG(hObject,handles)
    set(handles.cargandoPlot,'visible','off')
end



% --- Executes on button press in abrirnuevafigura.
function abrirnuevafigura_Callback(hObject, eventdata, handles)
% hObject    handle to abrirnuevafigura (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in botonExportar.
function botonExportar_Callback(hObject, eventdata, handles)
% hObject    handle to botonExportar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function rangotemporal_Callback(hObject, eventdata, handles)
MaxElegible = handles.MaxElegible;
Valor = str2double(get(hObject,'String'));
set(hObject,'String',num2str(floor(Valor)))
if str2double(get(hObject,'String')) > MaxElegible
    set(hObject,'String',num2str(floor(MaxElegible)))
end
if str2double(get(hObject,'String'))<2
    set(hObject,'String','2')
end


% --- Executes during object creation, after setting all properties.
function rangotemporal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangotemporal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in generarvideoboton.
function generarvideoboton_Callback(hObject, eventdata, handles)
% hObject    handle to generarvideoboton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of generarvideoboton



function limvideoms_Callback(hObject, eventdata, handles)
MaxElegible = handles.MaxElegible;
Valor = str2double(get(hObject,'String'));
set(hObject,'String',num2str(floor(Valor)))
if str2double(get(hObject,'String')) > MaxElegible
    set(hObject,'String',num2str(floor(MaxElegible)))
end
if str2double(get(hObject,'String'))==0
    set(hObject,'String','1')
end


% --- Executes during object creation, after setting all properties.
function limvideoms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to limvideoms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menuinicio_Callback(hObject, eventdata, handles)
% hObject    handle to menuinicio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function botonleerarchivos_Callback(hObject, eventdata, handles)
% Solicita los archivos
[NombreArchivo,Ruta] = uigetfile(...
    '*.*',...
    'Elige los archivos ''Post File''',...
    'MultiSelect', 'on');
if isempty(NombreArchivo); return; end
% Genera la ruta hacia los archivos
for i=1:size(NombreArchivo,2)
    handles.ListaArchivos(i) = strcat(Ruta,NombreArchivo(i));
end

UltimaCarpeta = find(Ruta==filesep,2,'last');
TituloCarpeta = Ruta(UltimaCarpeta(1):end);
handles.figure1.Name = ['EASEPostFile2Matlab',' (...',TituloCarpeta,')'];

% Importa los archivos
importarPostFiles(hObject,handles)

% Carga las nuevas variables
handles = guidata(hObject);

% Carga de datos inicial
% Extraer cada ID de fuente diferente. Sirve para tener identificada cada
% fuente y para conocer el numero de fuentes.
handles.Fuentes = unique(handles.SPLm(:,2));

% Indice de inicio de cada fuente
[~, handles.FuentesIni] = unique(handles.SPLm(:,2),'first');
% Indice de final de cada fuente
[~, handles.FuentesFin] = unique(handles.SPLm(:,2),'last');

% Indice de los datos inicial
handles.Index = 1;

% Array de nombres de receptores
handles.Receptores = handles.SPLm(...
    handles.FuentesIni(handles.Index):handles.FuentesFin(handles.Index),...
    3);

% Actualiza los valores mostrados en la aplicación
% Sala
set(handles.volumentexto,'String',num2str(handles.SPLm{handles.Index,15}))
set(handles.superficietexto,'String',num2str(handles.SPLm{handles.Index,16}))
set(handles.temperaturatexto,'String',num2str(handles.SPLm{handles.Index,17}))
set(handles.humedadtexto,'String',num2str(handles.SPLm{handles.Index,18}))
% Fuente
set(handles.posicionfuentetexto,'String',...
    [num2str(handles.SPLm{handles.Index,4}(1)),'x',...
    num2str(handles.SPLm{handles.Index,4}(2)),'x',...
    num2str(handles.SPLm{handles.Index,4}(3))])
set(handles.retardofuentetexto,'String',num2str(handles.SPLm{handles.Index,7}))
% Receptor
set(handles.posicionreceptortexto,'String',...
    [num2str(handles.SPLm{handles.Index,5}(1)),'x',...
    num2str(handles.SPLm{handles.Index,5}(2)),'x',...
    num2str(handles.SPLm{handles.Index,5}(3))])
set(handles.tiempollegadatexto,'String',num2str(handles.SPLm{handles.Index,9}))
set(handles.distanciafuentetexto,'String',num2str(handles.SPLm{handles.Index,6}))
set(handles.rayosrecibidostexto,'String',num2str(handles.SPLm{handles.Index,8}))
% Lista receptores
set(handles.listareceptores,'string',handles.Receptores)
% Lista fuentes
set(handles.listafuentes,'string',handles.Fuentes)

% Actualiza el valor maximo elegible
set(handles.ValMaxText,'String',num2str(floor(handles.MaxElegible)));


% Mostrar detalles del programa
set(handles.panelparametrossala,'visible','on')
set(handles.panelfuente,'visible','on')
set(handles.panellistafuentes,'visible','on')
set(handles.panellistareceptores,'visible','on')
set(handles.panelparametrosreceptor,'visible','on')
set(handles.paneltemporal,'visible','on')
set(handles.botonplotear,'visible','on')

% Actualizar variables en la aplicacion
guidata(hObject, handles);


% --- Executes on button press in botonAbsorcion.
function botonAbsorcion_Callback(hObject, eventdata, handles)
representaAbsorcion(handles)


% --- Executes on button press in teoriacorregida.
function teoriacorregida_Callback(hObject, eventdata, handles)
% hObject    handle to teoriacorregida (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of teoriacorregida



function Qdir_Callback(hObject, eventdata, handles)
%Valor = str2double(get(hObject,'String'));
%set(hObject,'String',num2str(floor(Valor)))
if str2double(get(hObject,'String')) < 1
    set(hObject,'String','1')
end


% --- Executes during object creation, after setting all properties.
function Qdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Qdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numIntentos_Callback(hObject, eventdata, handles)
Valor = str2double(get(hObject,'String'));
set(hObject,'String',num2str(floor(Valor)))
if str2double(get(hObject,'String')) < 10
    set(hObject,'String','10')
end


% --- Executes during object creation, after setting all properties.
function numIntentos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numIntentos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ecogramaBoton.
function ecogramaBoton_Callback(hObject, eventdata, handles)
% hObject    handle to ecogramaBoton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ecogramaBoton


% --------------------------------------------------------------------
function botonimportarmat_Callback(hObject, eventdata, handles)
% Solicita los archivos
[NombreArchivo,Ruta] = uigetfile(...
    '*.mat',...
    'Selecciona el archivo .mat');

Archivo = strcat(Ruta,NombreArchivo);

DatosCargados = load(Archivo);

UltimaCarpeta = find(Ruta==filesep,2,'last');
TituloCarpeta = Ruta(UltimaCarpeta(1):end);
handles.figure1.Name = ['EASEPostFile2Matlab',' (...',TituloCarpeta,NombreArchivo,')'];

handles.SPLm = DatosCargados.SPLm;
handles.MaxElegible = DatosCargados.MaxElegible;

% Carga de datos inicial
% Extraer cada ID de fuente diferente. Sirve para tener identificada cada
% fuente y para conocer el numero de fuentes.
handles.Fuentes = unique(handles.SPLm(:,2));

% Indice de inicio de cada fuente
[~, handles.FuentesIni] = unique(handles.SPLm(:,2),'first');
% Indice de final de cada fuente
[~, handles.FuentesFin] = unique(handles.SPLm(:,2),'last');

% Indice de los datos inicial
handles.Index = 1;

% Array de nombres de receptores
handles.Receptores = handles.SPLm(...
    handles.FuentesIni(handles.Index):handles.FuentesFin(handles.Index),...
    3);

% Actualiza los valores mostrados en la aplicación
% Sala
set(handles.volumentexto,'String',num2str(handles.SPLm{handles.Index,15}))
set(handles.superficietexto,'String',num2str(handles.SPLm{handles.Index,16}))
set(handles.temperaturatexto,'String',num2str(handles.SPLm{handles.Index,17}))
set(handles.humedadtexto,'String',num2str(handles.SPLm{handles.Index,18}))
% Fuente
set(handles.posicionfuentetexto,'String',...
    [num2str(handles.SPLm{handles.Index,4}(1)),'x',...
    num2str(handles.SPLm{handles.Index,4}(2)),'x',...
    num2str(handles.SPLm{handles.Index,4}(3))])
set(handles.retardofuentetexto,'String',num2str(handles.SPLm{handles.Index,7}))
% Receptor
set(handles.posicionreceptortexto,'String',...
    [num2str(handles.SPLm{handles.Index,5}(1)),'x',...
    num2str(handles.SPLm{handles.Index,5}(2)),'x',...
    num2str(handles.SPLm{handles.Index,5}(3))])
set(handles.tiempollegadatexto,'String',num2str(handles.SPLm{handles.Index,9}))
set(handles.distanciafuentetexto,'String',num2str(handles.SPLm{handles.Index,6}))
set(handles.rayosrecibidostexto,'String',num2str(handles.SPLm{handles.Index,8}))
% Lista receptores
set(handles.listareceptores,'string',handles.Receptores)
% Lista fuentes
set(handles.listafuentes,'string',handles.Fuentes)

% Actualiza el valor maximo elegible
set(handles.ValMaxText,'String',num2str(floor(handles.MaxElegible)));

% Mostrar detalles del programa
set(handles.panelparametrossala,'visible','on')
set(handles.panelfuente,'visible','on')
set(handles.panellistafuentes,'visible','on')
set(handles.panellistareceptores,'visible','on')
set(handles.panelparametrosreceptor,'visible','on')
set(handles.paneltemporal,'visible','on')
set(handles.botonplotear,'visible','on')

% Actualizar variables en la aplicacion
guidata(hObject, handles);


% --------------------------------------------------------------------
function exportarmenu_Callback(hObject, eventdata, handles)
% hObject    handle to exportarmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function exportmatboton_Callback(hObject, eventdata, handles)
[file,path,~] = uiputfile({...
    '*.mat','MAT-files (*.mat)';                             % 1 - JPEG
    },'Guardar datos','Datos.mat');

if file==0; return; end
SPLm=handles.SPLm;
MaxElegible = handles.MaxElegible;
% Guarda los datos importados de los TXT en .mat
save([path,file],'SPLm','MaxElegible')

% Si se ha calculado la teoria revisada se exportan el resto de datos
if isfield(handles,'Elate')
    
    % Coeficientes
    eL = handles.Elate;
    CL = handles.Clate;
    eE = handles.Eearly;
    CE = handles.Cearly;
    CD = handles.Cdirect;
    
    % Parametros
    SPL = handles.SPL; % Nivel de presion a 1 metro
    W = handles.W; % Potencia
    Q = handles.Q; % Directividad
    t = handles.t; % Rango temporal en segundos
    T = handles.T; % Tiempo de reverberación de Eyring
    c = handles.c; % Velocidad de sonido
    S = handles.S; % Superficie
    V = handles.V; % Volumen
    alpha = handles.alpha; % Coeficiente de absorcion medio
    Z = handles.Z; % Impedancia del aire
    rho = handles.rho; % Densidad del aire
    A = handles.A; % Absorcion equivalente de Eyring
    RTmid = mean([handles.SPLm{1,10}(8),...
        handles.SPLm{1,10}(11),...
        handles.SPLm{1,10}(14)]);
    
    % Curvas
    Dist = handles.Distancia; % Vector de distancia
    EaseU = handles.EaseU; % Curva del campo util (EASE)
    EaseP = handles.EaseP; % Curva del campo perjudicial (EASE)
    TeoU = handles.TeoU; % Curva del campo util (Teorico)
    TeoP = handles.TeoP; % Curva del campo perjudicial (Teorico)
    EaseD = handles.EaseD; % Curva del campo directo (EASE)
    EaseE = handles.EaseE; % Curva del campo early (EASE)
    EaseL = handles.EaseL; % Curva del campo late (EASE)
    TeoD = handles.TeoD; % Curva del campo directo (Teorico)
    TeoE = handles.TeoE; % Curva del campo early (Teorico)
    TeoL = handles.TeoL; % Curva del campo late (Teorico)
    
    if isfield(handles,'DistCorte')
        DistCorte = handles.DistCorte;
    else
        DistCorte = NaN;
    end
    if isfield(handles,'DistCorteTeo')
        DistCorteTeo = handles.DistCorteTeo;
    else
        DistCorteTeo = NaN;
    end
    
    save([path,file],...
        'eL','CL','eE','CE','CD',... % Coeficientes
        'SPL','W','Q','t','T','c','S','V','alpha','Z','rho','A','RTmid',... % Parametros
        'Dist','EaseU','EaseP','TeoU','TeoP',... % Curvas
        'EaseD','EaseE','EaseL','TeoD','TeoE','TeoL',... % Curvas 2
        'DistCorte','DistCorteTeo',... % Puntos de corte
        '-append')   
end
msgbox('Los datos se han exportado correctamente','Datos guardados','help','modal')


% --- Executes on button press in frec125.
function frec125_Callback(hObject, eventdata, handles)
% hObject    handle to frec125 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frec125


% --- Executes on button press in frec250.
function frec250_Callback(hObject, eventdata, handles)
% hObject    handle to frec250 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frec250


% --- Executes on button press in frec500.
function frec500_Callback(hObject, eventdata, handles)
% hObject    handle to frec500 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frec500


% --- Executes on button press in frec1000.
function frec1000_Callback(hObject, eventdata, handles)
% hObject    handle to frec1000 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frec1000


% --- Executes on button press in frec2000.
function frec2000_Callback(hObject, eventdata, handles)
% hObject    handle to frec2000 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frec2000


% --- Executes on button press in frec4000.
function frec4000_Callback(hObject, eventdata, handles)
% hObject    handle to frec4000 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frec4000


% --- Executes on button press in frec8000.
function frec8000_Callback(hObject, eventdata, handles)
% hObject    handle to frec8000 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frec8000


% --- Executes on button press in frecFull.
function frecFull_Callback(hObject, eventdata, handles)
% hObject    handle to frecFull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frecFull


% --- Executes on button press in botonPerjudicial.
function botonPerjudicial_Callback(hObject, eventdata, handles)
% hObject    handle to botonPerjudicial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of botonPerjudicial


% --- Executes on button press in botonEarly.
function botonEarly_Callback(hObject, eventdata, handles)
% hObject    handle to botonEarly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of botonEarly


% --- Executes on button press in botonDirecto.
function botonDirecto_Callback(hObject, eventdata, handles)
% hObject    handle to botonDirecto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of botonDirecto


% --- Executes on button press in trescampos.
function trescampos_Callback(hObject, eventdata, handles)
% Inhabilita botones que no tienen calculo con los tres campos separados
set(handles.claridadboton,'enable','off')
set(handles.definicionboton,'enable','off')
set(handles.sonoridadboton,'enable','off')
set(handles.ecogramaBoton,'enable','off')


% --- Executes on button press in doscampos.
function doscampos_Callback(hObject, eventdata, handles)
% Habilita los botones
set(handles.claridadboton,'enable','on')
set(handles.definicionboton,'enable','on')
set(handles.sonoridadboton,'enable','on')
set(handles.ecogramaBoton,'enable','on')
