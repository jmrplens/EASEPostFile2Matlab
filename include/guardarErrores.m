
function guardarErrores(Ucoef,Pcoef,Distancia,Dist_P,Dist_U)

a = confint(Ucoef);
b = confint(Pcoef);

% Resolucion de errores
x = linspace(min(Distancia),max(Distancia),8);

% Puntos de error del campo util
U = Ucoef.a*x.^(Ucoef.b);
Utop = a(2)*x.^(a(4))-U;
Udow = U-a(1)*x.^(a(3));
TUerrores = table(x',U',Utop',Udow','VariableNames',{'x','y','ymas','ymenos'});

writetable(TUerrores,'UtilErrores.dat','Delimiter','\t')

% Puntos campo perjudicial
% Puntos de error del campo util
P = Pcoef.p1*x+Pcoef.p2;
Ptop = (b(2)*x+(b(4)))-P;
Pdow = P-(b(1)*x+(b(3)));
TPerrores = table(x',P',Ptop',Pdow','VariableNames',{'x','y','ymas','ymenos'});

writetable(TPerrores,'PerjudicialErrores.dat','Delimiter','\t')
end