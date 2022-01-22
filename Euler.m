%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Programa - Método de Euler                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
% Engenharia Biomédica e Biofísica - Modelação e Simulação em Medicina
% Faculdade de Ciências Da Universidade de Lisboa, 2019
% Afonso Pourfarzaneh, 48249
% Ângelo Nunes, 48254
% Maria Lopes, 49941

clear all;
clear vars;
clc;

parametros;
%% Passo 0 - Atingir o equilíbrio do mecanismo NF-kB, sem RANKL (PR = 0)

% Definir intervalo de tempo para a integração da função NFK
start = 0;
fim = 54000;
delta_t = 10;                                                              % Passo temporal
n_pts = fim/delta_t;                                                       % Número de pontos
vetor_tempo = linspace(0,54000,n_pts);

% Definir matrizes para a derivada da variável y e para a variável em si
% que serão preenchidas pelo método de Euler - cada coluna corresponde a
% uma equação diferencial.

y = zeros(n_pts,15);                                                      
dy = zeros(n_pts,15);                                                      
PR = 0;                                                                    % Ausência de RANKL
y(1,:) = y0;                                                               % Condições iniciais de cada equação diferencial

% Método de Euler
for i = 1:n_pts-1
    dy(i,:) = NFK(vetor_tempo,y(i,:),PR);
    y(i+1,:) = y(i,:)+delta_t.*dy(i,:);
end

%% Passo 1  - Progresso do mecanismo NF-kB, ainda sem RANKL (PR = 0) 

% Definir novo intervalo de tempo para a continuação da função NFK
start1 = 0;
fim1 = 18000;
delta_t1 = 10;
n_pts1 = fim1/delta_t1;
vetor_tempo1 = linspace(start1,fim1,n_pts1);

% Definir novas matrizes

dy1 = zeros(n_pts1,15);
y1 = zeros(n_pts1,15);
x1 = size(y);
y1(1,:) = y(x1(1),:);                                                      % Atualizar novas condições iniciais que correspondem aos últimos 15 valores obtidos no passo anterior

% Método de Euler
for i = 1:n_pts1-1
    dy1(i,:) = NFK(vetor_tempo1,y1(i,:),PR);
    y1(i+1,:) = y1(i,:)+delta_t1.*dy1(i,:);
end
%% Passo 2 - Introdução da Remodelação Óssea com adição de 4 equações diferenciais e presença de RANKL (PR = 1)

% Definir intervalo de tempo para a integração da função NFK_bone
start2 =18000;                                                             %começa onde a NFK acaba
fim2 = 36000;
delta_t2 = 1;
n_pts2 = (fim2 - start2)/delta_t2 ;

% Repartição do intervalo de tempo para implementar as mudanças rápidas e
% lentas das EDOs (discutido no relatório)
start2_1 = 36000;
fim2_1 = 136000;
delta_t2_1 = 100;
n_pts2_1 = (fim2_1 - start2_1)/delta_t2_1;

% vetor tempo total
vetor_tempo2 = linspace(start2,fim2,n_pts2);
vetor_tempo2_1 = linspace(start2_1,fim2_1,n_pts2_1);
vetor_tempo2 = [vetor_tempo2 , vetor_tempo2_1];
n_pts2 = n_pts2 + n_pts2_1;

% Definir matrizes, agora com 19 colunas devido á adição das 4 equações
% diferenciais da remodelação óssea
dy2 = zeros(n_pts2,19);
y2 = zeros(n_pts2,19);
x2 = size(y1);
PR = 1;
y2(1,:) = [y1(x2(1),:),7.626e-4,6.332e-4,1.049e-4,100];                    % Atualizar e adicionar condições iniciais para as 19 EDOs

% Método de Euler
for i = 1:n_pts2-1
    dy2(i,:) = NFK_bone(vetor_tempo2,y2(i,:),PR);
    y2(i+1,:) = y2(i,:)+delta_t2.*dy2(i,:);
end
%% Passo 3 - Começo da simulação do efeito da injeção/inibição da proteína RANKL

% Definir intervalo de tempo 1, em que se integra a função NFK_bone
start3 =0;
fim3 = 18000;
delta_t3 = 1;
n_pts3 = (fim3 - start3)/delta_t3 ;
vetor_tempo3 = linspace(start3,fim3,n_pts3);  

% Definir matrizes
dy3 = zeros(n_pts3,19);
y3 = zeros(n_pts3,19);
x3 = size(y2);
y3(1,:) = y2(x3(1),:);                                                     %Condições iniciais são os últimos valores obtidos no passo 2

% Método de Euler
for i = 1:n_pts3-1
    dy3(i,:) = NFK_bone(vetor_tempo3,y3(i,:),PR);
    y3(i+1,:) = y3(i,:)+delta_t3.*dy3(i,:);
end
%% Passo 4 - Injeção/Inibição da RANKL usando a função NFK_bone_RANKL

% Continuação do intervalo de tempo definido anteriormente
start4 =18000;
fim4 = 64000;
delta_t4 = 1;
n_pts4 = (fim4 - start4)/delta_t4 ;
vetor_tempo4 = linspace(start4,fim4,n_pts4);

% Definir matrizes
dy4 = zeros(n_pts4,19);
y4 = zeros(n_pts4,19);
x4 = size(y3);
y4(1,:) = y3(x4(1),:);                                                     % Condições iniciais

% Pedir ao utilizador se deseja simular a injeção ou inibição
decisao = input(" Inibição (0) ou injeção(1): ");

if decisao == 1
    for i = 1:n_pts4-1
        inibicao = 1;                                                      % Inibição de 10%
        injecao = 20;
        dy4(i,:) = NFK_bone_RANKL(vetor_tempo4,y4(i,:),PR,injecao,inibicao);
        y4(i+1,:) = y4(i,:)+delta_t4.*dy4(i,:);
    end
elseif decisao == 0
    for i = 1:n_pts4-1
        injecao = 0;                                                       % Injeção de 10 mM de RANKL
        inibicao = 0.9;
        dy4(i,:) = NFK_bone_RANKL(vetor_tempo4,y4(i,:),PR,injecao,inibicao);
        y4(i+1,:) = y4(i,:)+delta_t4.*dy4(i,:);
    end
end

%% Passo 5 - Fim da simulação 

% Finalizar intervalo de tempo para a simulação
start5 = 64000;
fim5 = 100000;
delta_t5 = 1;
n_pts5 = (fim5 - start5)/delta_t5 ;
vetor_tempo5 = linspace(start5,fim5,n_pts5);

% Definir matrizes
dy5 = zeros(n_pts5,19);
y5 = zeros(n_pts5,19);
x5 = size(y4);
y5(1,:) = y4(x5(1),:);

% Método de Euler
for i = 1:n_pts5-1
    dy5(i,:) = NFK_bone(vetor_tempo5,y5(i,:),PR);
    y5(i+1,:) = y5(i,:)+delta_t5.*dy5(i,:);
end
%% Gráficos


% Junção dos resultados e dos intervalos de tempo
Y = [y3;y4;y5];
T = [vetor_tempo3';vetor_tempo4';vetor_tempo5'];
T = T/3600;                                                                % conversão do tempo de segundos para horas

% Para cada gráfico os valores encontram-se normalizados, ou seja os
% valores no eixo das ordenadas são adimensionais 

figure(1)
subplot(1,3,1);plot(T,Y(:,7)./Y(1,7));grid on;title('NFkB');
subplot(1,3,2);plot(T,Y(:,15)./Y(1,15));grid on;title('TFs');
subplot(1,3,3);
plot(T,Y(:,16)./Y(1,16));hold on;
plot(T,Y(:,17)./Y(1,17),'r-');hold on;
plot(T,Y(:,18)./Y(1,18),'k-.');title('OBp, OBa, OCa');
legend('OBp','OBa','OCa');grid on;

figure(2)
subplot(1,2,1);plot(T,(Y(:,17)./Y(1,17))./(Y(:,18)./Y(1,18)),'r-');
title('Rácio OBa:OCa');
grid on; 
subplot(1,2,2);plot(T,Y(:,19)./Y(1,19)*100);grid on;title('BV');

