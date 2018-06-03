function representarG(hObject,handles)
% Reprsentacion del parametro sonoridad

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

% Temperatura
temp = handles.SPLm{handles.Index,17};
% Velocidad del sonido
c = 331.4 + 0.6*temp;
% Densidad del aire. Presion Atm / (constante gas*Temp
% absoluta)
rho = 101325 / (287.058*(273.15+temp));
% Impedancia acustica del aire
Z = rho*c;
% Nivel de presion total de la fuente a 1 metro
SPL = 10*log10(sum(10.^(handles.SPLm{handles.Index,13}/10)));
% Nivel de potencia de la fuente (SPL a 1 metro). Se supone Q=1
SWL = SPL - 10*log10(1/(4*pi*1));
% Potencia de la fuente en vatios
W = 10.^(SWL/10)*10^-12;

%% Representacion grafica del campo Util y el Perjudicial

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


% Claridad EASE
% Muestra los puntos
plot(Dist,10*log10( (10.^(Urec/10) + 10.^(Prec/10)) ./ 10^((Dcoef.a*10^Dcoef.b)/10) ) ,'^b','MarkerFaceColor','b'),hold on
% Muestra la curva
Curv(1)=plot(Distplot,10*log10( (10.^(Uplot/10) + 10.^(Pplot/10)) ./ 10^((Dcoef.a*10^Dcoef.b)/10) ),'b','LineWidth',1.2);

% Añade las rejillas de las graficas
grid on
grid minor

% Leyenda
leyenda{1} = 'G - EASE';

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
    
    Curv(2) = plot(Distplot,10*log10( (10.^(Pdeplot/10) + 10.^(Plateplot/10)) ./ (Fac * ( W./(4*pi*10.^2)*fitiD.Cd)) ),'-.r','LineWidth',1);
    
    % Leyendas
    leyenda{2} = 'G - Teórico';
    
end

lgdw=legend(Curv,leyenda);
lgdw.FontSize = 11;
%lgdw.Location = 'northeastoutside';
% Si se utiliza Matlab 2014b o anterior no añade el titulo a la leyenda ya
% que no es compatible
if ~verLessThan('matlab','8.7')
    title(lgdw,'Curvas')
end
ylabel('G (dB)')
xlabel('Distancia (m)')
title({'G';...
    ['\color{blue} ','Fuente: ',sprintf(': %s',Fuente(1))]});

hold off

