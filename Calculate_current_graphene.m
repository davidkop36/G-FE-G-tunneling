function I1=Calculate_current_graphene(kt1,kt2,k_vi,Q_mis)
flag = 2;

if nargin == 3

for ii=1:size(kt1,1)
    for jj=1:size(kt1,2)
        k = linspace(kt1(ii,jj),kt2(ii,jj),1000);

        I1(ii,jj) = trapz(k, abs(k-k_vi(ii,jj)).*abs(k).*(1+k.^2 + (k - k_vi(ii,jj) ).^2  ) ./ ( (1 + k_vi(ii,jj)^2)^(3/2) .* (1 + (2*k-k_vi(ii,jj)).^2 ).^(3/2)    ));


    end
end


else
    if flag==1
    phi1 = linspace(0,2*pi,100)';
    phi2 = linspace(0,2*pi,100)';
    unity_ma = ones(size(phi1));
    
    for ii=1:size(kt1,1)
        for jj=1:size(kt1,2)
            k = linspace(kt1(ii,jj),kt2(ii,jj),100)';
            %full_exp = ((kt2(ii,jj)>kt1(ii,jj))*2-1)  *abs(k).*abs(k-k_vi(ii,jj)  ) ./ (1 + k.^2 + (k-k_vi(ii,jj)).^2 + Q_mis.^2 + 2 *tensorprod(abs(k).*abs(k-k_vi(ii,jj)),(tensorprod(cos(phi1),cos(phi2)) + tensorprod(sin(phi1),sin(phi2) ) )) + abs(Q_mis)*(tensorprod(abs(k), tensorprod(cos(phi1),unity_ma)  ) +  tensorprod(abs(k-k_vi(ii,jj)), tensorprod(unity_ma,cos(phi2))  )   )    ).^2;
            full_exp = abs(k).*abs(k-k_vi(ii,jj)  ) ./ (1 + k.^2 + (k-k_vi(ii,jj)).^2 + Q_mis.^2 + 2 *tensorprod(abs(k).*abs(k-k_vi(ii,jj)),(tensorprod(cos(phi1),cos(phi2)) + tensorprod(sin(phi1),sin(phi2) ) )) +2* abs(Q_mis)*(tensorprod(abs(k), tensorprod(cos(phi1),unity_ma)  ) +  tensorprod(abs(k-k_vi(ii,jj)), tensorprod(unity_ma,cos(phi2))  )   )    ).^2;
            full_exp = squeeze(full_exp);
            I1(ii,jj) = trapz(k,trapz(phi1,trapz(phi2, full_exp,3),2),1);
    
        end
    end
    
    
    
    else
            phi = linspace(0,2*pi,105);
            unity_ma = ones(size(phi));
                for ii=1:size(kt1,1)
                    for jj=1:size(kt1,2)

                        k = linspace(kt1(ii,jj),kt2(ii,jj),100)';
                        a1 = 1 + k.^2 + (k-k_vi(ii,jj)).^2 + Q_mis.^2;
                        c1 = 2*abs(Q_mis.*k);
                        b1 = 2*abs(k.*(k-k_vi(ii,jj)));
                        d1 = 2*abs(Q_mis.*(k-k_vi(ii,jj)));
                        full_exp = 2*pi*abs(k).*abs(k-k_vi(ii,jj)  ).*abs(a1+c1*cos(phi)) ./ (abs(a1.^2+c1.^2-b1.^2-d1.^2 + 2*(a1.*c1-b1.*d1)*cos(phi) )).^(3/2) ;
                        I1(ii,jj) = trapz(k,trapz(phi, full_exp,2),1);
                
                    end
                end

    end

end
