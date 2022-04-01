function [delta_xi,EMF]  = zhiqianChen(g1,g2,omega_e,R_s,Ld,Lq,xi,u_alphabeta,i_alphabeta)
    delta_xi = [0;0];
    EMF = [0;0];
    %p¦Î=    A_0 ¦Î   +   A_1 i   +   A_2 u
    %A_0=   G/L_d+¦Ø_re  J 
    %A_1=   G¦Ø_re  J+G(R_s-?L¦Ø_re J+G)/L_d
    %A_2=   - G/L_d
    G = [g1 -g2;
        g2 g1];
    I = [1 0;
        0 1];
    J = [0 -1;
        1 0];
    
    A_0 = G/Ld + omega_e*J;
    A_1 = G*omega_e*J + G*(R_s-(Ld-Lq)*omega_e*J +G)/Ld;
    A_2 = -G/Ld;
    delta_xi = A_0*xi + A_1*i_alphabeta + A_2*u_alphabeta;
    EMF = xi + G*i_alphabeta;

end
