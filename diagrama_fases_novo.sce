clc
clear

function psat = antoine(T)
    A = [4.02832 ; 5.24677];
    B = [1268.636 ; 1598.673];
    C = [-56.199 ; -46.424];
    //psat = 10^(A - B./(T+C)); //psat em bar T em Kelvin
    A1 = [6.89677 ; 8.1122];
    B1 = [1264.9  ; 1592.864];
    C1 = [216.544 ; 226.184]; 
    psat = 10^(A1 - B1./(T-273.15 + C1))*0.0013604323308265306
endfunction

function coef_ativ = gama_rk(x,T)
    denset = 0.7956;
    maset = 40;
    denshep = 0.68867;
    mashep = 100;
    R = 1.987;
    v1=mashep/denshep; 
    v2=maset/denset;
    v2 = 58
    V = [v1 ; v2]
    x2 = 1 - x; 
    x = [x ; x2];
    azinho = [617 ; 2096];
    //azinho = [935.36;1593.97];
    Gama = [(V(2)/V(1))*exp(-azinho(1)/(1.987*T));(V(1)/V(2))*exp(-azinho(2)/(1.987*T))];
    MiE_RT1 = -log(x(1) + Gama(1)*x(2)) + x(2)*((Gama(1)/(x(1) + Gama(1)*x(2)))-(Gama(2)/(x(2)+Gama(2)*x(1))));
    MiE_RT2 = -log(x(2) + Gama(2)*x(1)) - x(1)*((Gama(1)/(x(1) + Gama(1)*x(2)))-(Gama(2)/(x(2)+Gama(2)*x(1))));
    MiE_RT = [MiE_RT1;MiE_RT2];
    //coef_ativ = exp(MiE_RT)
    //coef_ativ = [1;1]
    teta = azinho;
    x1 = x(1);
    x2 = x(2);
    lambda_12 = v2/v1*exp(-teta(1)/(R*(T))) //aqui, T deve estar em Kelvin
    lambda_21 = v1/v2*exp(-teta(2)/(R*(T)))
    beta_aux = (lambda_12/(x1+lambda_12*x2) - lambda_21/(lambda_21*x1 + x2))
    gama1 = exp(-log(x1+lambda_12*x2)+x2*beta_aux)
    gama2 = exp(-log(x2+lambda_21*x1)-x1*beta_aux);
    coef_ativ = [gama1 ; gama2];
    //coef_ativ = [1;1]
endfunction


function [gama1,gama2] = wilson(T,x1,x2)
    denset =  0.7956
    maset = 40
    denshep = 0.68867
    mashep = 100  
    R = 1.987
    v1=mashep/denshep; 
    //v2=maset/denset;
    v2 = 58;
    teta = [746 ; 658];
    lambda_12 = v2/v1*exp(-teta(1)/(R*(T))) //aqui, T deve estar em Kelvin
    lambda_21 = v1/v2*exp(-teta(2)/(R*(T)))
    beta_aux = (lambda_12/(x1+lambda_12*x2) - lambda_21/(lambda_21*x1 + x2))
    gama1 = exp(-log(x1+lambda_12*x2)+x2*beta_aux)
    gama2 = exp(-log(x2+lambda_21*x1)-x1*beta_aux)
endfunction


function F = sistema_bolha(P,x,teta)
    y = teta(1);
    T = teta(2);
    coef_ativ = gama_rk(x,T);
    psat = antoine(T);
    F1 = x*coef_ativ(1)*psat(1) - y*P;
    F2 = (1 - x)*coef_ativ(2)*psat(2) - (1-y)*P;
    F = [F1 ; F2];
endfunction

function F = sistema_orvalho(P,y,teta)
    x = teta(1);
    T = teta(2);
    coef_ativ = gama_rk(x,T);
    psat = antoine(T);
    F1 = x*coef_ativ(1)*psat(1) - y*P;
    F2 = (1 - x)*coef_ativ(2)*psat(2) - (1-y)*P;
    F = [F1 ; F2];
endfunction

function F = sistema_azeo(P,teta)
    x = teta(1);
    T = teta(2);
    coef_ativ = gama_rk(x,T);
    psat = antoine(T);
    F1 = coef_ativ(1)*psat(1) - P;
    F2 = coef_ativ(2)*psat(2) - P;
    F = [F1 ; F2];
endfunction

function F = sistema_banco(P,teta,q)
    x = teta(1);
    T = teta(2);
    coef_ativ = gama_rk(x,T);
    psat = antoine(T);
    F1 = coef_ativ(1)*psat(1) - P - q(1);
    F2 = coef_ativ(2)*psat(2) - P - q(2);
    F = [F1 ; F2];
endfunction

function J = jacobiana_bolha(P,x,teta)
    J = [];
    h = 1e-5;
    for k = 1:2
        teta_adv = teta
        teta_adv(k) = teta_adv(k) + h;
        der = (sistema_bolha(P,x,teta_adv) - sistema_bolha(P,x,teta))/h;
        J = [J der];
    end
endfunction

function J = jacobiana_orvalho(P,y,teta)
    J = [];
    h = 1e-5;
    for k = 1:2
        teta_adv = teta
        teta_adv(k) = teta_adv(k) + h;
        der = (sistema_orvalho(P,y,teta_adv) - sistema_orvalho(P,y,teta))/h;
        J = [J der];
    end
endfunction

function J = jacobiana_azeo(P,teta)
    J = [];
    h = 1e-5;
    for k = 1:2
        teta_adv = teta
        teta_adv(k) = teta_adv(k) + h;
        der = (sistema_azeo(P,teta_adv) - sistema_azeo(P,teta))/h;
        J = [J der];
    end
endfunction

function J = jacobiana_banco(P,teta,q)
    J = [];
    h = 1e-5;
    for k = 1:2
        teta_adv = teta
        teta_adv(k) = teta_adv(k) + h;
        der = (sistema_banco(P,teta_adv,q) - sistema_banco(P,teta,q))/h;
        J = [J der];
    end
endfunction

function teta_bolha = newton_bolha(P,x,teta0)
    erro = 1
    while erro > 1e-8
        J = jacobiana_bolha(P,x,teta0);
        F = sistema_bolha(P,x,teta0);
        novoteta = teta0 - 0.5*inv(J)*F
        erro = norm(novoteta - teta0);
        teta0 = novoteta
    end
    teta_bolha = novoteta;
endfunction


function teta_orvalho = newton_orvalho(P,y,teta0)
    erro = 1
    while erro > 1e-8
        J = jacobiana_orvalho(P,y,teta0);
        F = sistema_orvalho(P,y,teta0);
        novoteta = teta0 - 0.5*inv(J)*F
        erro = norm(novoteta - teta0);
        teta0 = novoteta
    end
    teta_orvalho = novoteta;
endfunction

function teta_azeo = newton_azeo(P,teta0)
    erro = 1
    while erro > 1e-8
        J = jacobiana_azeo(P,teta0);
        F = sistema_azeo(P,teta0);
        novoteta = teta0 - inv(J)*F
        erro = norm(novoteta - teta0);
        teta0 = novoteta
    end
    teta_azeo = novoteta;
endfunction

function teta_azeo = newton_banco(P,teta0,q)
    erro = 1
    while erro > 1e-8
        J = jacobiana_banco(P,teta0,q);
        F = sistema_banco(P,teta0,q);
        novoteta = teta0 - inv(J)*F
        erro = norm(novoteta - teta0);
        teta0 = novoteta
    end
    teta_azeo = novoteta;
endfunction

teta_azeo = newton_azeo(1.0133,[.3; 65+274])

pause

// gerando as curvas de ponto de bolha

vx_b = [];
vT_b = [];
P = 0.8;
teta0 = [.1;450];

for x = 0.01:0.001:0.999
    teta_bolha = newton_bolha(P,x,teta0);
    vx_b = [vx_b ; x];
    vT_b = [vT_b ; teta_bolha(2)];
    teta0 = teta_bolha;
end

plot(vx_b,vT_b-273.15,'b-')
xtitle('Wilson Model','$x_1, y_1$','$T(ºC)$') ;
vx_o = [];
vT_o = [];
teta0 = [.1;350];

for y = 0.01:0.001:0.999
    teta_orvalho = newton_orvalho(P,y,teta0);
    vx_o = [vx_o ; y];
    vT_o = [vT_o ; teta_orvalho(2)];
    teta0 = teta_orvalho;
end

plot(vx_o,vT_o-273.15,'r-')

experimento(:,1)=[0.0172742  
    0.0263820  
    //0.0272851  
    0.0293095  
    //0.0355040  
    //0.0423977  
    0.0448166  
    0.0573863  
    0.1002261  
    //0.1249792  
    0.1666204  
    0.2028871  
    0.2251639  
    0.2827184  
    //0.2946279  
    //0.3504181
    0.3706239  
    0.4551522
    0.7977160
    //0.8056999
    0.9053357  
    //0.9256985  
    //0.9279876
    0.9279876  
    //0.9582398  
    //0.9653529
    //0.9845673  
    //0.9869946  
    //0.9894275
    ]

experimento(:,2)=[0.0256672  
    0.0638794  
    //0.0957119  
    0.0864437  
    //0.1081857  
    //0.1404002  
    0.1712420  
    0.2028871  
    0.2682215  
    //0.2773022  
    0.3214778  
    0.3180105  
    0.3188743  
    0.3180105  
    //0.2780705  
    //0.2586680
    0.2842822  
    0.2811618
    0.3302897
    //0.2922124
    0.4621636  
    //0.5151488  
    //0.8036969
    0.4780102  
    //0.9511765  
    //0.8680246
    //0.6190385  
    //0.6462330  
    //0.7761423 
    ];

experimento(:,3)=273.15 + [72.4
70.125
//69.25
69.52
//69
//68.09
67.13
65.99
65.1
//64.72
64.49
64.34
64.31
64.25
//64.19
//64.05
64.15
63.83
64.41
//65.52
69.19
//74.88
//85.02
70.33
//90.34
//88.37
//78.02
//79.59
//84.09
]

plot(experimento(:,1),experimento(:,3)-273.15,'bo')

plot(experimento(:,2),experimento(:,3)-273.15,'bx')
