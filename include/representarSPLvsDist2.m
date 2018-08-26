function representarSPLvsDist2(hObject,handles)

axes(handles.plotfig)

Rango = handles.Rango;
SPLmTot = handles.SPLm(:,19); % Matriz con los rayos sin procesar
SPLm = handles.SPLm(:,20); % Matriz con valores cada ms y por octavas
Distancia = cell2mat(handles.SPLm(:,6))'; % Vector de distancias fuente-receptor
Distplot = min(Distancia):0.01:max(Distancia); % Vector de distancias para las curvas
Fuente = string(handles.SPLm(:,2)');
Receptor = string(handles.SPLm(:,3)');


% Comprobar que bandas de frecuencias se han elegido
Bandas = [];
fsel = [];
if handles.frec125.Value==1
    Bandas = [Bandas,1];
    fsel = [fsel,125];
end
if handles.frec250.Value==1
    Bandas = [Bandas,2];
    fsel = [fsel,250];
end
if handles.frec500.Value==1
    Bandas = [Bandas,3];
    fsel = [fsel,500];
end
if handles.frec1000.Value==1
    Bandas = [Bandas,4];
    fsel = [fsel,1000];
end
if handles.frec2000.Value==1
    Bandas = [Bandas,5];
    fsel = [fsel,2000];
end
if handles.frec4000.Value==1
    Bandas = [Bandas,6];
    fsel = [fsel,4000];
end
if handles.frec8000.Value==1
    Bandas = [Bandas,7];
    fsel = [fsel,8000];
end
if isempty(Bandas); Bandas = 1:7; end
if isempty(fsel); fsel = 10^3*(2.^((-3:3))); end

%% Procesado de datos
% Separacion e integracion de los niveles antes y despues del rango
Dneto = zeros(1,numel(SPLm));
Eneto = zeros(1,numel(SPLm));
Lneto = zeros(1,numel(SPLm));
for i=1:numel(SPLm)
    % Direct - Early - Late EXTRACTION
    Dbruto = SPLm{i}(1,Bandas); % Niveles de 0 a 1 ms
    Ebruto = SPLm{i}(2:Rango,Bandas); % Niveles de 1 ms a tiempo de integracion
    Lbruto = SPLm{i}(Rango+1:end,Bandas); % Niveles de tiempo de integracion a infinito
    
    % Ponderacion A
    %     f = fsel; % Frecuencias seleccionadas
    %     ponderacionA = 20*log10(...
    %         12194^2*f.^4 ./...
    %         ((f.^2+20.6^2).*sqrt((f.^2+107.7^2).*(f.^2+737.9^2)).*(f.^2+12194^2)))+2;
    %     Dbruto = Dbruto + ponderacionA;
    %     Ebruto = Ebruto + ponderacionA;
    %     Lbruto = Lbruto + ponderacionA;
    
    % Suma de niveles para los dos rangos
    Dneto(i) = 10*log10(sum(10.^(Dbruto(:)/10)));
    Eneto(i) = 10*log10(sum(10.^(Ebruto(:)/10)));
    Lneto(i) = 10*log10(sum(10.^(Lbruto(:)/10)));
end

% Campo util es la suma del campo directo y el early
Uneto = 10*log10(10.^(Dneto/10)+10.^(Eneto/10));

%% Representacion grafica del campo Util y el Perjudicial
colores = colormap(jet(21));
% Indices de color
iCU = 1;
iCUt = 3;
iCL = 18;
iCLt = 21;
%% Puntos y curva de 0ms a Valor
% Ordenar valores de distancia de menor a mayor manteniendo su valor de SPl
% asociado
Dist_U = sortrows([Distancia',Uneto']);
% Si tiene distancias duplicadas/iguales se promedia el nivel de las dos
% distancias iguales y elimina el duplicado
[C,~,idx] = unique(Dist_U(:,1),'stable');
val = accumarray(idx,Dist_U(:,2),[],@(x) 10*log10(mean(10.^(x/10))));
Dist_U = [C,val];
% Añadir por extrapolacion valor en la posicion de 0 metros hasta la posicion del
% receptor mas cercano
%V0_0toVal=interp1(Dist_U(:,1),Dist_U(:,2),0,'linear','extrap');
%Dist_U = [[0;V0_0toVal]';Dist_U];
% Separar los valores para graficar y convertirlo en vector fila
Dist = Dist_U(:,1);
Urec = Dist_U(:,2);
% Crear curva potencial. Obtiene coeficientes (Ucoef) y datos de ajuste
% (Ures)
[Ucoef,Ures,Uoutput]=fit(Dist,Urec,'power1');
Uplot=Ucoef.a*Distplot.^Ucoef.b;
% Muestra los puntos
plot(Dist,Urec,'s','Color',colores(iCU,:),'MarkerFaceColor',colores(iCU,:))
hold on
% Muestra los errores
Uerr = confint(Ucoef);
Ucurvarec = Ucoef.a*Dist.^Ucoef.b;
errorbar(Dist,Ucoef.a*Dist.^Ucoef.b,... % Curva
    (Uerr(2)*Dist.^Uerr(4)) - Ucurvarec,... % Error positivo
    Ucurvarec - (Uerr(1)*Dist.^Uerr(3)),... % Error negativo
    'b','LineStyle','none')
% Muestra la curva
Curv(1)=plot(Distplot,Uplot,'Color',colores(iCU,:),'LineWidth',1.2);

%% Puntos y curva de valor a infinito ms
% Ordenar valores de distancia de menor a mayor manteniendo su valor de SPl
% asociado
Dist_P = sortrows([Distancia',Lneto']);
% Si tiene distancias duplicadas/iguales se promedia el nivel de las dos
% distancias iguales y elimina el duplicado
[C,~,idx] = unique(Dist_P(:,1),'stable');
val = accumarray(idx,Dist_P(:,2),[],@(x) 10*log10(mean(10.^(x/10))));
Dist_P = [C,val];
% Añadir por extrapolacion valor en la posicion de 0 metros hasta la posicion del
% receptor mas cercano
%V0_0toVal=interp1(Dist_P(:,1),Dist_P(:,2),0,'linear','extrap');
%Dist_P = [[0;V0_0toVal]';Dist_P];
% Separar los valores para graficar y convertirlo en vector fila
Dist = Dist_P(:,1);
Prec = Dist_P(:,2);
% Obtener coeificentes curva polinomica de grado 1
[Pcoef,Pres,Poutput] = fit(Dist,Prec,'poly1');
% Crea la curva
Pplot = Pcoef.p1*Distplot+Pcoef.p2;
% Muestra los puntos
plot(Dist,Prec,'^','Color',colores(iCL,:),'MarkerFaceColor',colores(iCL,:))
% Muestra los errores
Perr = confint(Pcoef);
Pcurvarec = Pcoef.p1*Dist+Pcoef.p2;
errorbar(Dist,Pcoef.p1*Dist+Pcoef.p2,... % Curva
    (Perr(2)*Dist+Perr(4)) - Pcurvarec,... % Error positivo
    Pcurvarec - (Perr(1)*Dist+Perr(3)),... % Error negativo
    'r','LineStyle','none')
% Muestra la curva
Curv(2)=plot(Distplot,Pplot,'Color',colores(iCL,:),'LineWidth',1.2);

% Añade las rejillas de las graficas
grid on
grid minor

% Obtener punto de cruce
CortesInd = find(abs(Uplot-Pplot)<=(0.01));
if ~isnan(CortesInd)
    DistCorte = Distplot(CortesInd(1));
    text(DistCorte,Uplot(CortesInd(1))+1,...
        {'EASE',strcat('\bf\color{black}',[sprintf('%4.3f',DistCorte),' m'])},...
        'VerticalAlignment','bottom','HorizontalAlignment','left')
end

%% Guarda un .dat para generar graficas en el LaTex
%guardarErrores(Ucoef,Pcoef,Distancia,Dist_P,Dist_U)
%%

% Obtener numero de ceros decimales de las pendientes, sumarle 1 y convertir a
% string, para mostrar tantos numeros decimales en las funciones que se ven en
% la leyenda
P11 = num2str(fix(abs(log10(abs(Ucoef.b))))+2);
P21 = num2str(fix(abs(log10(abs(Pcoef.p1))))+2);

% Colores en formato texto para el color de las funciones
ColU = [num2str(colores(iCU,1)),',',num2str(colores(iCU,2)),',',num2str(colores(iCU,3))];
ColL = [num2str(colores(iCL,1)),',',num2str(colores(iCL,2)),',',num2str(colores(iCL,3))];

% Leyenda
leyenda{1} = sprintf(['\\bf0 %s %d ms\\rm   R^2_{adj} = %4.2f\n\\color[rgb]{%s}y = %4.2f·x^{%4.',P11,'f} \n \\color{white}.'],'a',Rango,Ures.adjrsquare,ColU,Ucoef.a,Ucoef.b);
leyenda{2} = sprintf(['\\bf%d ms %s\\rm   R^2_{adj} = %4.2f\n\\color[rgb]{%s}y = %4.',P21,'f·x+%4.2f \n \\color{white}.'],Rango,'a infinito',Pres.adjrsquare,ColL,Pcoef.p1,Pcoef.p2);

%% Teoria corregida
% Si se ha marcado el calculo de la teoria revisada corregida se calcula y
% se añade a la grafica
if get(handles.teoriacorregida,'value')==1
    %%% PARAMETROS
    % Nivel de presion total de la fuente a 1 metro
    SPL = 10*log10(sum(10.^(handles.SPLm{handles.Index,13}/10)));
    % Volumen (m3)
    V = handles.SPLm{handles.Index,15};
    % Superficie (m2)
    S = handles.SPLm{handles.Index,16};
    % Mean Free Path
    mfp = 4*V/S;
    % Absorcion media
    alpha = mean(handles.SPLm{handles.Index,11});
    % Absorción equivalente de Eyring
    A = -S*log(1-alpha);
    % Temperatura de la sala
    temp = handles.SPLm{handles.Index,17};
    % Velocidad del sonido
    c = 331.4 + 0.6*temp;
    % Densidad del aire. Presion Atm / (constante gas*Temp
    % absoluta)
    rho = 101325 / (287.058*(273.15+temp));
    % Impedancia acustica del aire
    Z = rho*c;
    % Presion de referencia
    pref = 2*10^(-5);
    % Factor de conversion
    Fac = Z/(pref^2);
    % Directividad de la fuente
    Q = str2double(get(handles.Qdir,'string'));
    % Nivel de potencia de la fuente (SPL a 1 metro). Se supone Q=1
    SWL = SPL - 10*log10(1/(4*pi*1));
    % Potencia de la fuente en vatios
    W = 10.^(SWL/10)*10^-12;
    % Tiempo de reverberacion de Eyring
    T = 6*log(10)/c * mfp * 1/(-log(1-alpha));
    % Rango de integracion temporal (s)
    t = handles.Rango * 10^(-3);
    
    %% Obtencion de la curvas directo y early
    %% Curva del campo directo
    % Ordenar valores de distancia de menor a mayor manteniendo su valor de SPl
    % asociado
    Dist_D = sortrows([Distancia',Dneto']);
    % Si tiene distancias duplicadas/iguales se promedia el nivel de las dos
    % distancias iguales y elimina el duplicado
    [C,~,idx] = unique(Dist_D(:,1),'stable');
    val = accumarray(idx,Dist_D(:,2),[],@(x) 10*log10(mean(10.^(x/10))));
    Dist_D = [C,val];
    % Añadir por extrapolacion valor en la posicion de 0 metros hasta la posicion del
    % receptor mas cercano
    %V0_0toVal=interp1(Dist_U(:,1),Dist_U(:,2),0,'linear','extrap');
    %Dist_U = [[0;V0_0toVal]';Dist_U];
    % Separar los valores para graficar y convertirlo en vector fila
    Dist = Dist_D(:,1);
    Drec = Dist_D(:,2);
    % Crear curva potencial. Obtiene coeficientes (Ucoef) y datos de ajuste
    % (Ures)
    [Dcoef,Dres]=fit(Dist,Drec,'power1');
    Dplot=Dcoef.a*Distplot.^Dcoef.b;
    
    %% Curva del campo early
    % Ordenar valores de distancia de menor a mayor manteniendo su valor de SPl
    % asociado
    Dist_E = sortrows([Distancia',Eneto']);
    % Si tiene distancias duplicadas/iguales se promedia el nivel de las dos
    % distancias iguales y elimina el duplicado
    [C,~,idx] = unique(Dist_E(:,1),'stable');
    val = accumarray(idx,Dist_E(:,2),[],@(x) 10*log10(mean(10.^(x/10))));
    Dist_E = [C,val];
    % Añadir por extrapolacion valor en la posicion de 0 metros hasta la posicion del
    % receptor mas cercano
    %V0_0toVal=interp1(Dist_U(:,1),Dist_U(:,2),0,'linear','extrap');
    %Dist_U = [[0;V0_0toVal]';Dist_U];
    % Separar los valores para graficar y convertirlo en vector fila
    Dist = Dist_E(:,1);
    Erec = Dist_E(:,2);
    % Crear curva potencial. Obtiene coeficientes (Ucoef) y datos de ajuste
    % (Ures)
    [Ecoef,Eres]=fit(Dist,Erec,'power1');
    Eplot=Ecoef.a*Distplot.^Ecoef.b;
    
    % Para no mostrar errores cuando se repite el calculo de coeficientes
    warning ('off','all');
    
    %% Calculo de coeficiente para el campo directo
    for ii = 1:str2double(get(handles.numIntentos,'string')) % 100 intentos de obtener los coeficientes
        try % Intenta el calculo con puntos de inicio random
            opD = fitoptions('Method','NonlinearLeastSquares',...
                'Robust','LAR',...
                'Lower',[0,0],...   % Minimo valor de los coeficientes
                'Upper',[20,20],... % Maximo valor de los coeficientes
                'MaxFunEvals',3000,...
                'MaxIter',3000,...
                'TolFun',10^-3);
            fD = fittype('Fac * W*Q./(4*pi*dist.^2)*Cd',...
                'dependent',{'y'},...
                'independent',{'dist'},...
                'problem', {'Q','W','Fac'},...
                'coefficients',{'Cd'},...
                'options',opD);
            [fitiD,gofD,outputD] = fit(Distplot',10.^(Dplot'/10),fD,'problem',{Q,W,Fac});
            break;
        catch % Si no se ha podido obtener unos coeficientes ejecuta lo siguiente
            if ii==str2double(get(handles.numIntentos,'string'))
                hold off
                ed = errordlg('No se ha podido obtener coeficientes (Directo), intentalo de nuevo','Error');
                set(ed, 'WindowStyle', 'modal');
                uiwait(ed);
                return;
            end
        end
    end
    
    
    %% Calculo de coeficientes epsilon y Cl para campo late
    for ii = 1:str2double(get(handles.numIntentos,'string')) % 100 intentos de obtener los coeficientes
        try % Intenta el calculo con puntos de inicio random
            fL = fittype('Fac * (4*W)/(S*(-log(1-alpha))) * exp(-(13.82*(dist/c+t)*el/T)) * Cl',...
                'dependent',{'y'},...
                'independent',{'dist'},...
                'problem', {'W','S','alpha','t','T','c','Fac'},...
                'coefficients',{'el','Cl'});
            [fitiL,gofL,outputL] = fit(Distplot',10.^(Pplot'/10),fL,'problem',{W,S,alpha,t,T,c,Fac});
            break;
        catch % Si no se ha podido obtener unos coeficientes ejecuta lo siguiente
            if ii==str2double(get(handles.numIntentos,'string'))
                hold off
                ed = errordlg('No se ha podido obtener coeficientes (Late), intentalo de nuevo','Error');
                set(ed, 'WindowStyle', 'modal');
                uiwait(ed);
                return;
            end
        end
    end
    
    %% Calculo de coeficientes epsilon y Ce para campo early
    for ii = 1:str2double(get(handles.numIntentos,'string')) % 100 intentos de obtener los coeficientes
        try % Intenta el calculo con puntos de inicio random
            opE = fitoptions('Method','NonlinearLeastSquares',...
                'Robust','LAR',...
                'MaxFunEvals',3000,...
                'MaxIter',3000,...
                'TolFun',10^-3);
            fE = fittype('Fac * ((4*W)./(S*(-log(1-alpha))*dist) .* (exp(-(13.82*(dist/c)*ee/T))*Ce - exp(-(13.82*(dist/c+t)*el/T)) * Cl))',...
                'dependent',{'y'},...
                'independent',{'dist'},...
                'problem', {'W','S','alpha','T','t','c','Fac','el','Cl'},...
                'coefficients',{'ee','Ce'},...
                'options',opE);
            [fitiE,gofE,outputE] = fit(Distplot',10.^(Eplot'/10),fE,'problem',{W,S,alpha,T,t,c,Fac,fitiL.el,fitiL.Cl});
            break;
        catch % Si no se ha podido obtener unos coeficientes ejecuta lo siguiente
            if ii==str2double(get(handles.numIntentos,'string'))
                hold off
                ed = errordlg('No se ha podido obtener coeficientes (Early), intentalo de nuevo','Error');
                set(ed, 'WindowStyle', 'modal');
                uiwait(ed);
                return;
            end
        end
    end
    warning ('on','all');
    
    PDplot = 10*log10(Fac * ( W*Q./(4*pi*Distplot.^2)*fitiD.Cd));
    Peplot = 10*log10(Fac *(...
        (4*W)./(S*(-log(1-alpha))*Distplot) .* ...
        (exp(-(13.82*(Distplot/c)*fitiE.ee/T))*fitiE.Ce - exp(-(13.82*(Distplot/c+t)*fitiL.el/T))*fitiL.Cl)));
    Pdeplot = 10*log10(10.^(PDplot/10)+10.^(Peplot/10));
    Plateplot = 10*log10(Fac .* (4*W)./(S*(-log(1-alpha))) .* exp(-(13.82.*(Distplot./c+t).*fitiL.el/T)) .* fitiL.Cl);
    
    Curv(3) = plot(Distplot,Pdeplot,'-.','Color',colores(iCUt,:),'LineWidth',1);
    Curv(4) = plot(Distplot,Plateplot,'-.','Color',colores(iCLt,:),'LineWidth',1);
    
    % Colores en formato texto para el color de las funciones
    ColUteo = [num2str(colores(iCUt,1)),',',num2str(colores(iCUt,2)),',',num2str(colores(iCUt,3))];
    ColLteo = [num2str(colores(iCLt,1)),',',num2str(colores(iCLt,2)),',',num2str(colores(iCLt,3))];
    
    % Leyendas
    leyenda{3} = sprintf('\\bf0 %s %d ms\\rm   Teoria revisada corregida\n \\color[rgb]{%s}\\epsilon_E = %4.3f | C_E = %4.3f | R^2_{adj} = %4.2f \n C_D = %4.3f | R^2_{adj} = %4.2f \n \\color{white}.','a',Rango,ColUteo,fitiE.ee,fitiE.Ce,gofE.adjrsquare,fitiD.Cd,gofD.adjrsquare);
    leyenda{4} = sprintf('\\bf%d ms %s\\rm   Teoria revisada corregida\n \\color[rgb]{%s}\\epsilon_L = %4.3f | C_L = %4.3f | R^2_{adj} = %4.2f',Rango,'a infinito',ColLteo,fitiL.el,fitiL.Cl,gofL.adjrsquare);
    
    
    % Obtener punto de cruce
    CortesInd = find(abs(Pdeplot-Plateplot)<=(0.01));
    if ~isnan(CortesInd)
        DistCorteTeo = Distplot(CortesInd(1));
        text(DistCorteTeo,Pdeplot(CortesInd(1))-1,...
            {'Teórica',strcat('\bf\color{black}',[sprintf('%4.3f',DistCorteTeo),' m'])},...
            'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
    
    % Almacenar variables de interes para despues exportarlas dentro del
    % .mat
    
    % Coeficientes
    handles.Elate = fitiL.el;
    handles.Clate = fitiL.Cl;
    handles.Eearly = fitiE.ee;
    handles.Cearly = fitiE.Ce;
    handles.Cdirect = fitiD.Cd;
    
    % Parametros
    handles.SPL = SPL; % Nivel de presion a 1 metro
    handles.W = W; % Potencia
    handles.Q = Q; % Directividad
    handles.t = t; % Rango temporal en segundos
    handles.T = T; % Tiempo de reverberación de Eyring
    handles.c = c; % Velocidad de sonido
    handles.S = S; % Superficie
    handles.V = V; % Volumen
    handles.alpha = alpha; % Coeficiente de absorcion medio
    handles.Z = Z; % Impedancia del aire
    handles.rho = rho; % Densidad del aire
    handles.A = A; % Absorcion equivalente de Eyring
    
    % Curvas
    handles.Distancia = Distplot; % Vector de distancia
    handles.EaseU = Uplot; % Curva del campo util (EASE)
    handles.EaseP = Pplot; % Curva del campo perjudicial (EASE)
    handles.TeoU = Pdeplot; % Curva del campo util (Teorico)
    handles.TeoP = Plateplot; % Curva del campo perjudicial (Teorico)
    handles.EaseD = Dplot; % Curva del campo directo (EASE)
    handles.EaseE = Eplot; % Curva del campo early (EASE)
    handles.EaseL = Pplot; % Curva del campo late (EASE)
    handles.TeoD = PDplot; % Curva del campo directo (Teorico)
    handles.TeoE = Peplot; % Curva del campo early (Teorico)
    handles.TeoL = Plateplot; % Curva del campo late (Teorico)
    
    if exist('DistCorte','var')
        handles.DistCorte = DistCorte;
    end
    if exist('DistCorteTeo','var')
        handles.DistCorteTeo = DistCorteTeo;
    end
    % Actualizar variables en la aplicacion
    guidata(hObject, handles);
    
end

lgdw=legend(Curv,leyenda);
lgdw.FontSize = 11;
%lgdw.Location = 'northeastoutside';
% Si se utiliza Matlab 2014b o anterior no añade el titulo a la leyenda ya
% que no es compatible
if ~verLessThan('matlab','8.7')
    title(lgdw,'Curvas')
end
ylabel('Nivel de presión acústica (dB)')
xlabel('Distancia (m)')
title({['SPL según rango de tiempo',sprintf(' - %d ms',Rango)];...
    ['\color{blue} ','Fuente: ',sprintf(': %s',Fuente(1))]});

hold off

% Asegurar el tamaño y posicion del plot
set(gca, 'units','normalized', 'outerposition',[-0.073 -0.027 0.753 1.028]);
set(gca, 'units','normalized', 'position',[0.037 0.077 0.562 0.854]);

