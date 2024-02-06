function [dYdt] = Ca_model(~,Y, pars)
% Cao et al. (2013) Ca model

c=Y(1);
B=Y(2);

Jrelease=pars(1);
No=pars(2);
Jleak=pars(3);
BT=pars(4);
Vs=pars(5);
Ks=pars(6);
k_on=pars(7);
k_off=pars(8);

Jserca=Vs*c/(c+Ks);

dYdt=[Jrelease*No+Jleak-Jserca+k_off*(BT-B)-k_on*c*B,...
        k_off*(BT-B)-k_on*c*B
    ];
end

