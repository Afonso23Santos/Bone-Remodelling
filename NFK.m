function dy_dt = NFK(~,vetor_y,PR)
parametros;

dy_dt = zeros(15,1);                                                       %Definir vetor de outputs
 
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

end

