function [Area]     = get_IntArea_guess(vd,vc,peak_abs)

Phi_D_vo    = (2/vd)*sqrt(log(2)/pi);   
a           = sqrt(log(2))*vc/vd;      

V_vo        = exp(a^2)*erfc(a);       
Voigt_vo    = Phi_D_vo*V_vo;         

Area        = peak_abs/Voigt_vo;
