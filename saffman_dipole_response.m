function ux = saffman_dipole_response(xs,ys,zs,rs,bx,by,bz,params)
% NOTE WE STILL NEED TO HAVE THE SAME SET OF ARGUMENTS IN ORDER FOR THIS TO
% WORK EVEN THOUGH zs, rs, and bz are unused

assert(params.Gelastic==params.G2D,'Gelastic should be equal to G2D to run saffman')

F = params.F;
d = params.d;
G2D = params.G2D;
Gint = params.Gint;
saff_length = G2D/(Gint); 


xs1 = xs - bx*d/2;
ys1 = ys - by*d/2;
%zs1 = zs - bz*d/2;
rs1 = sqrt(xs1.^2+ys1.^2);

xs2 = xs + bx*d/2;
ys2 = ys + by*d/2;
%zs2 = zs + bz*d/2;
rs2 = sqrt(xs2.^2+ys2.^2);

kr1 = rs1./saff_length;
kr2 = rs2./saff_length;
        
f1 = StruveH0(kr1)-StruveH1(kr1)./kr1-1/2*(bessely(0,kr1)-bessely(2,kr1))+2./(pi*kr1.^2);
f2 = StruveH0(kr2)-StruveH1(kr2)./kr2-1/2*(bessely(0,kr2)-bessely(2,kr2))+2./(pi*kr2.^2);
        
g1 = StruveH0(kr1) - 2*StruveH1(kr1)./kr1 + bessely(2,kr1) + 4./(pi*kr1.^2);
g2 = StruveH0(kr2) - 2*StruveH1(kr2)./kr2 + bessely(2,kr2) + 4./(pi*kr2.^2);
        
ux= (F.*(f1.*bx - g1.*((xs1.*xs1./rs1.^2).*bx + (xs1.*ys1./rs1.^2).*by)) ... 
            -F.*(f2.*bx - g2.*((xs2.*xs2./rs2.^2).*bx + (xs2.*ys2./rs2.^2).*by)))/4; % overall factor of 4
                                                                                     % note there should be an overall
                                                                                     % factor of 1/G_2D
                                                                                     % -- this is implemented by choosing G_0 = G_2D

        %ux = (params.F*params.d./(8*pi*rs.^2)).*(3*(bx.*xs./rs + by.*ys./rs + bz.*zs./rs).^2-1).*(xs./rs);