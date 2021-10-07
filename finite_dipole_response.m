function ux = finite_dipole_response(xs,ys,zs,rs,bx,by,bz,params)

F = params.F;
d = params.d;

xs1 = xs - bx*d/2;
ys1 = ys - by*d/2;
zs1 = zs - bz*d/2;
rs1 = sqrt(xs1.^2+ys1.^2+zs1.^2);

xs2 = xs + bx*d/2;
ys2 = ys + by*d/2;
zs2 = zs + bz*d/2;
rs2 = sqrt(xs2.^2+ys2.^2+zs2.^2);


ux = (((rs1.^(-1)).*(F*bx + (xs1.*xs1./rs1.^2).*(F.*bx) + (xs1.*ys1./rs1.^2).*F.*by + (xs1.*zs1./rs1.^2).*F.*bz) ...
    -(rs2.^(-1)).*(F*bx + (xs2.*xs2./rs2.^2).*(F.*bx) + (xs2.*ys2./rs2.^2).*(F.*by) + (xs2.*zs2./rs2.^2).*(F.*bz))))/(8*pi);
%ux = (params.F*params.d./(8*pi*rs.^2)).*(3*(bx.*xs./rs + by.*ys./rs + bz.*zs./rs).^2-1).*(xs./rs);