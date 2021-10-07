function ux = ideal_dipole_response(xs,ys,zs,rs,bx,by,bz,params)
ux = (params.F*params.d./(8*pi*rs.^2)).*(3*(bx.*xs./rs + by.*ys./rs + bz.*zs./rs).^2-1).*(xs./rs);