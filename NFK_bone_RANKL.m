function dy_dt = NFK_bone_RANKL(~,vetor_y,PR,injecao,inibicao)

parametros;
%% TGF

Alfa=1;                                                                    %TGF-Beta armazenado na matriz óssea
D_TGF=2.00e+2;                                                             %Taxa de degradação de TGF-Beta.
S_TGF=0;                                                                   %Dose externa de TGF-Beta.
TGF=(Alfa*Kr*vetor_y(18)+S_TGF)/D_TGF;                                     %Concentração de TGF.
%% Ligação do TGF
K_D1_TGF=4.8250e-04;                                                       %Coeficiente de Ativação relativo á ligação TGF-Beta nos OBu
K_D2_TGF=2.1900e-04;                                                       %Coeficiente de Repressão relativo á ligação TGF-Beta nos OBp
K_D3_TGF=9.3304e-05;                                                       %Coeficiente de Ativação relativo á ligação TGF-Beta nos OCa

%Funções hill

PI_TGF_act_OBu=TGF/(K_D1_TGF+TGF);                                         %Função hill da estimulação do TGF-Beta sobre os OBu
PI_TGF_rep_OBp=1/(1+TGF/K_D2_TGF);                                         %Função hill da repressão do TGF-Beta sobre os OBp
PI_TGF_act_OCa=TGF/(K_D3_TGF+TGF);                                         %Função hill da estimulação do TGF-Beta sobre os OCa
%% PTH
Beta_PTH=9.74e+2;                                                          %Taxa de sintetização de PTH sistémico
P_PTH_d=0;                                                                 %Dose externa de PTH.
D_PTH=3.84e+2;                                                             %Taxa de degradação de PTH.
PTH=(Beta_PTH+P_PTH_d)/D_PTH;                                              %Concentração de PTH

%% OPG
K_D2_PTH = 2.21e-1;                                                        %Coef. de Repressão para a produção de OPG relacionado com ligação do RANKL
Beta_2_OPG = 5.42e+6;                                                      %Taxa de produção mínima de OPG por OBa
P_OPG_d = 0;                                                               %Dose externa de OPG
OPG_max = 7.98e+2;                                                         %Maxima concentração de OPG possível
D_OPG = 4.16;                                                              %Taxa de degradação de OPG

%Função Hill

PI_PTH_rep_OB_a = 1/(1+PTH/K_D2_PTH);                                      %Função Hill da repressão do PTH sobre os OBa

Numerador = Beta_2_OPG*vetor_y(17)*PI_PTH_rep_OB_a+P_OPG_d;
Denominador = (Beta_2_OPG*vetor_y(17)*PI_PTH_rep_OB_a)/OPG_max+D_OPG;
OPG = Numerador/Denominador;

%% IL6
K_M_TGF_IL6_act = 2.9e-3;                                                  %Coef. de ativação do IL6 pelo TGF
P_IL6_d = 0;                                                               %Dose externa de IL6
Beta_IL6 = 1.2e+7;                                                         %Taxa de produção mínima de IL6
IL6max = 8.04e-1;                                                          %Maxima concentração de IL6
D_IL6 = 4.99e+1;                                                           %Taxa de degradação de IL6

%Função Hill
PI_TGF_act_IL6 = TGF/(TGF+K_M_TGF_IL6_act);                                %Função Hill da ativação do TGF sobre o IL6
Part1 = P_IL6_d+Beta_IL6*OBu*PI_TGF_act_IL6;
Part2 = Beta_IL6*OBu*PI_TGF_act_IL6/IL6max+D_IL6;
IL6 = Part1/Part2;
%% RANKL
P_RANKL_d = 0 + injecao;                                                   %Dose externa de RANKL.
K_A_OPG = 5.68e-2;                                                         %Constante da taxa de ligação RANKL - OPG
K_A_RANK = 7.19e-2;                                                        %Constante da taxa de ligação RANKL - RANK
Beta_RANKL = 8.25e+5*inibicao;                                             %Produção endogena de RANKL por cada OBp
R_RANKL = 3e+6;                                                            %Maximo de RANKL expressado na superficie dos OBs precursores 
K_M_IL6_RANKL_act = 8.8;
K_D4_PTH = 2.09e+2;
RANK = 1.28e+1;
D_RANKL=4.16;                                                              %Taxa de degradação do RANKL

%funções hill
PI_IL6_act_RANKL=IL6/(IL6+K_M_IL6_RANKL_act);                              
PI_PTH_act_RANKL=PTH/(K_D4_PTH+PTH);
PI_ligands_RANKL=PI_IL6_act_RANKL+PI_PTH_act_RANKL- ...
PI_IL6_act_RANKL*PI_PTH_act_RANKL;

RANKL_Part1=P_RANKL_d+Beta_RANKL*vetor_y(16);
RANKL_Part2=1+K_A_OPG*OPG+K_A_RANK*RANK;
RANKL_Part3=Beta_RANKL/(R_RANKL*PI_ligands_RANKL)+D_RANKL;
RANKL=RANKL_Part1/(RANKL_Part2*RANKL_Part3);

%% Ligação do RANKL
K_D8_RANKL=4.12e+1;                                                        %Coef. de ativação relacionado com a ligação RANKL-RANK

%Função Hill
PI_RANKL_act_OCp=RANKL/(K_D8_RANKL+RANKL);

%% PR=0 RANKL off, PR=1 RANKL on

K_D9_RANKL=2.5;

%Função Hill
PI_RANKL=RANKL/(K_D9_RANKL+RANKL);
PR=PR*PI_RANKL; 
%% TFs

K_TFS=0.00065;

%Função Hill
PI_TFs_act_OCp=vetor_y(15)/(K_TFS+vetor_y(15));

%% Equações
dy_dt = zeros(19,1);                                                       %vetor outputs


dy_dt(1)=kprod-kdeg*vetor_y(1)-PR*k1*vetor_y(1);                           %IKK neutro   

dy_dt(2)=PR*k1*vetor_y(1)-k3*vetor_y(2)-PR*k2*vetor_y(2)*vetor_y(8) ...
-kdeg*vetor_y(2)-a2*vetor_y(2)*vetor_y(10)+t1*vetor_y(4)-a3*vetor_y(2) ...
*vetor_y(13)+t2*vetor_y(5);                                                %IKK ativo                                                                                    

dy_dt(3)=k3*vetor_y(2)+PR*k2*vetor_y(2)*vetor_y(8)-kdeg*vetor_y(3);        %IKK inativo  

dy_dt(4)=a2*vetor_y(2)*vetor_y(10)-t1*vetor_y(4);                          %concentração citoplasmática de (IKK|IkBa)  

dy_dt(5)=a3*vetor_y(2)*vetor_y(13)-t2*vetor_y(5);                          %concentração citoplasmática de (IKK|IkBa|NFkB) 

dy_dt(6)=c6a*vetor_y(13)-a1*vetor_y(6)*vetor_y(10)+ ...
t2*vetor_y(5)-i1*vetor_y(6);                                               %concentração citoplasmática de NFkB livre

dy_dt(7)=i1*kv*vetor_y(6)-a1*vetor_y(11)*vetor_y(7);                       %concentração nuclear de NFkB livre

dy_dt(8)=c4*vetor_y(9)-c5*vetor_y(8);                                      %concentração citoplasmática de A20

dy_dt(9)=c2+c1*vetor_y(7)-c3*vetor_y(9);                                   %transcrição de A20

dy_dt(10)=-a2*vetor_y(2)*vetor_y(10)-a1*vetor_y(10)*vetor_y(6)+ ...
c4a*vetor_y(12)-c5a*vetor_y(10)-i1a*vetor_y(10)+e1a*vetor_y(11);           %concentração citoplasmática de IkBa livre

dy_dt(11)=-a1*vetor_y(11)*vetor_y(7)+i1a*kv*vetor_y(10)- ...
e1a*kv*vetor_y(11);                                                        %concentração nuclear de IkBan livre

dy_dt(12)=c2a+c1a*vetor_y(7)-c3a*vetor_y(12);                              %transcrição de IkB

dy_dt(13)=a1*vetor_y(10)*vetor_y(6)-c6a*vetor_y(13)- ...
a3*vetor_y(2)*vetor_y(13)+e2a*vetor_y(14);                                 %concentração citoplasmática de(IkBa|NFkB) 

dy_dt(14)=a1*vetor_y(11)*vetor_y(7)-e2a*kv*vetor_y(14);                    %concentração nuclear de (IkBa|NFkB)

dy_dt(15)=c2c+c1c*vetor_y(7)-c3c*vetor_y(15);                              %concentração dos TFs - Fatores de transcrição regulados pelo NFKB
                                                          
dy_dt(16) = D_OBu*OBu*PI_TGF_act_OBu-D_OBp*PI_TGF_rep_OBp*vetor_y(16);     %Variação da concentração dos OBp

dy_dt(17) = D_OBp*PI_TGF_rep_OBp*vetor_y(16)-A_OBa*vetor_y(17);            %Variação da concentração dos OBa

dy_dt(18) = D_OCp*OCp*PI_TFs_act_OCp-A_OCa*vetor_y(18)*PI_TGF_act_OCa;     %Variação da concentração dos OCa

dy_dt(19) = -Kr*vetor_y(18)+Kf*vetor_y(17);                                %Variação da volume ósseo (BV)

end

