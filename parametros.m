%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Par�metros                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                           
%% Condi��es iniciais para o passo 0

y0=zeros(1,15);                                                            %Condi��es iniciais          
y0(13)=0.06;                                                               %n�vel inicial de concentra��o citoplasm�tica de(IkBa|NFkB)
%% Equa��es 1-15 -> mecanismo NFkB 

kv=5;                                                                      %r�cio volume citoplasm�tico/volume nuclear
  
% A20                         
                          
c1 = 5e-7;                                                                 
c2 = 0;                                                                    
c3 = 4e-4;                                                                  
c4 = 0.5;                                                                  
c5 = 3e-4;                                                                 
 
% IKK neutro e IKK ativo
 
k1=2.5e-3;             
k2=2.5e-5;                        
k3=1.5e-3;             
kprod=2.5e-5;          
kdeg=1.25e-4;          

% IkB alpha 
 
a1=0.5;                
a2=0.2;                
a3=1.;                 
t1=0.1;                  
t2=0.1;                

c1a=5e-7;              
c2a=0;                  
c3a=0.0004;            
c4a=0.005*100;          
c5a=0.0001;            
c6a=0.00002;           
 
i1=2.5e-3;             
e2a=0.01;              
i1a=1e-3;              
e1a=5e-4;              
  
 
% TFS
 
c1c=5e-7;              
c2c=0;                 
c3c=4e-4;              
 
%% Equa��es 16-19 -> Remodela��o �ssea
 
% Eq. 16 - Varia��o de OBp
D_OBu = 3.24e+2;                                                           %Taxa de diferencia��o de OBs progenitores
D_OBp = 3.67e-1;                                                           %Taxa de diferencia��o de OBs precursores
OBu = 3.27e-6;                                                             %Concentra��o de OBs n�o comprometidos (uncommited)

% Eq. 17 - Varia��o de OBa 
A_OBa = 3.00e-1;                                                           %Taxa de elimina��o de OBs ativos
 
% Eq. 18 - Varia��o de OCa
OCp = 1.28e-3;                                                             %Concentra��o de OCs precursores
D_OCp = 1.730e-1;                                                          %Taxa de diferencia��o de OCs precursores
A_OCa = 1.2;                                                               %Taxa de elimina��o de OCs ativos
 
% Eq. 19 - Varia��o de volume �sseo (BV) 
Kr=1.923690e+2;                                                            %Taxa relativa de reabsor��o �ssea
Kf=3.3148658888e+1;                                                        %Taxa relativa de forma��o �ssea