    function V = electrode(rho,I,xe,ye,x,y )
        q = rho * I;
        d = sqrt((xe-x).^2+(ye-y).^2);
        V = q./(4.*pi.*d);
    end