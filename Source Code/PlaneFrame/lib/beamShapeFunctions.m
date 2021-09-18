function [N,dNdxi] = beamShapeFunctions(xi)
    N     = [ ((xi - 1)^2*(xi + 2))/4
              ((xi - 1)^2*(xi + 1))/4
             -((xi + 1)^2*(xi - 2))/4
              ((xi - 1)*(xi + 1)^2)/4 ];
    dNdxi = [ (xi - 1)^2/4 + ((2*xi - 2)*(xi + 2))/4
              (xi - 1)^2/4 + ((2*xi - 2)*(xi + 1))/4
             -(xi + 1)^2/4 - ((2*xi + 2)*(xi - 2))/4
              (xi + 1)^2/4 + ((2*xi + 2)*(xi - 1))/4 ];
end