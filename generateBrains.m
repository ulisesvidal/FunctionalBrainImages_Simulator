function [stack_generated,stack_generated_masked,stack_labels,score_generated]=generateBrains(N,cdf,xmesh,SCORE,COEFF,mn,dim,mask,class)
% generateBrains() genera N cerebros de acuerdo a la clase especificada.
% 
% INPUTS:
%          N  - número de cerebros a generar
%       class - clase de cerebros a generar ('0 = normal','1= MCI','2= AD') 
%         cdf - matriz cuyas columnas son las cdfs generadas de cada score
%       xmesh - grids sobre los que se computan las cdfs

% OUTPUTS:
%     stack_generated - array 4D con N imágenes generadas de clase 'class' 


%% Generación de nuevos scores.
% Usamos la función emprand() para generar N puntos de acuerdo a la
% distribución estadística de cada score.

for s=1:size(SCORE,2)
    %Eliminamos muestras de las CDFs de cada score repetidas:
    auxcdf=cdf(:,s);
    rep= find(diff(auxcdf)==0);
    auxcdf(rep)=[];
    auxxmesh=xmesh(:,s);
    auxxmesh(rep)=[];

    %Scores generados:
    sc_generated_kde(:,s) = emprand(N,auxxmesh,auxcdf);
    
end

%% Reconstrucción de cerebros
%Montamos cerebros (estimación cdf: KDE):
X = sc_generated_kde*COEFF';
[m,~]=size(X);
stack_generated_masked= X+repmat(mn,m,1); %sumamos la media
stack_generated = zeros(N, dim(1), dim(2), dim(3),'uint8');
stack_generated(:,mask) =uint8(stack_generated_masked); 

%% Creamos vector de etiquetas ( NOR=0, MCI=1, AD=2)
stack_labels=[];
if class==0
    stack_labels=zeros(N,1);
elseif class==1
    stack_labels=zeros(N,1);
else
    stack_labels=ones(N,1)*2;
end

% clear cdf cdf0 cdf1 cdf2 xmesh xmesh0 xmesh1 xmesh2 class COEFF dim inf mask mn N SCORE stack_all inf
score_generated=sc_generated_kde;
end
