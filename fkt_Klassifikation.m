function [A, B1, B2, R] = fkt_Klassifikation(J, xWert, yWert)

    [VJ,DJ] = eig(J);
    lambda = DJ(1,1);
    mu = DJ(2,2);
    
    if mu == lambda
        if mu * lambda == 0
            [A, B1, B2, R] = fkt_entar_Knoten_EW0(INPUTS);
        else
            if lambda > 0
                if rank(VJ) == 1
                    [A, B1, B2, R] = fkt_entar_instab_Knoten_Ndiag(INPUTS);
                else
                    [A, B1, B2, R] = fkt_entar_instab_Knoten_diag(INPUTS);
                end
            else
                if rank(VJ) == 1
                    [A, B1, B2, R] = fkt_entar_stab_Knoten_Ndiag(INPUTS);
                else
                    [A, B1, B2, R] = fkt_entar_stab_Knoten_diag(INPUTS);
                end
            end
        end
    
        
        
    else %lambda ungleich mu
        if lambda * mu == 0
            if mu == 0
                if lambda > 0
                    [A, B1, B2, R] = fkt_lin_Expansion_mu0(INPUTS);
                else
                    [A, B1, B2, R] = fkt_lin_Kontraktion_mu0(INPUTS);
                end
            else
                if mu > 0
                    [A, B1, B2, R] = fkt_lin_Expansion_lambda0(INPUTS);
                else 
                    [A, B1, B2, R] = fkt_lin_Kontraktion_lambda0(INPUTS);
                end
            end
            
        else
            if abs(imag(lambda)) < 0.001
                if lambda * mu < 0
                    [A, B1, B2, R] = fkt_Sattel(INPUTS);
                else
                    if lambda > 0
                        [A, B1, B2] = fkt_Knoten(INPUTS);
                        R = -1; %instabil
                    else
                        [A, B1, B2] = fkt_Knoten(INPUTS);
                        R = 1; %stabil
                    end
                end
                
            else %lambda imaginaer
                if abs(real(lambda)) < 0.001
                    A = fkt_Zentrum(INPUTS);
                    R = 1;
                    B1 = 0;
                    B2 = 0;
                else
                    if real(lambda) > 0
                        A = fkt_Strudel(INPUTS);
                        R = -1; %instabil
                        B1 = 0;
                        B2 = 0;
                    else
                        A = fkt_Strudel(INPUTS);
                        R = 1; %stabil
                        B1 = 0;
                        B2 = 0;
                    end
                end
            end
        end
    end

end