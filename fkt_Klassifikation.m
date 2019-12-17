function [A, B1, B2] = fkt_Klassifikation(J, xWert, yWert)

    [VJ,DJ] = eig(J);
    lambda = DJ(1,1);
    mu = DJ(2,2);
    
    if mu == lambda
        if mu * lambda == 0
            [A, B1, B2] = fkt_entar_Knoten_EW0(INPUTS);
        else
            if lambda > 0
                if rank(VJ) == 1
                    [A, B1, B2] = fkt_entar_instab_Knoten_Ndiag(INPUTS);
                else
                    [A, B1, B2] = fkt_entar_instab_Knoten_diag(INPUTS);
                end
            else
                if rank(VJ) == 1
                    [A, B1, B2] = fkt_entar_stab_Knoten_Ndiag(INPUTS);
                else
                    [A, B1, B2] = fkt_entar_stab_Knoten_diag(INPUTS);
                end
            end
        end
    
        
        
    else %lambda ungleich mu
        if lambda * mu == 0
            if mu == 0
                if lambda > 0
                    [A, B1, B2] = fkt_lin_Expansion_mu0(INPUTS);
                else
                    [A, B1, B2] = fkt_lin_Kontraktion_mu0(INPUTS);
                end
            else
                if mu > 0
                    [A, B1, B2] = fkt_lin_Expansion_lambda0(INPUTS);
                else 
                    [A, B1, B2] = fkt_lin_Kontraktion_lambda0(INPUTS);
                end
            end
            
        else
            if abs(imag(lambda)) < 0.001
                if lambda * mu < 0
                    [A, B1, B2] = fkt_Sattel(INPUTS);
                else
                    if lambda > 0
                        [A, B1, B2] = fkt_instab_Knoten(INPUTS);
                    else
                        [A, B1, B2] = fkt_stab_Knoten(INPUTS);
                    end
                end
                
            else %lambda imaginaer
                if abs(real(lambda)) < 0.001
                    [A, B1, B2] = fkt_Zentrum(INPUTS);
                else
                    if real(lambda) > 0
                        [A, B1, B2] = fkt_instab_Strudel(INPUTS);
                    else
                        [A, B1, B2] = fkt_stab_Strudel(INPUTS);
                    end
                end
            end
        end
    end

end