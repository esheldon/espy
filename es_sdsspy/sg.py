def avg_gri(flux_g, ivar_g, 
            flux_r, ivar_r, 
            flux_i, ivar_i):
    """

    if you don't want to use an object in a particular and,  you should
    set the ivar to zero *before* calling this function

    fluxes are clipped between [0.001,1.e4], which is magnitude [12.5,30]
    """

    # clip the fluxes on the high and low end
    # this is mag between 12.5 and 30
    flux_g = flux_g.clip(0.001, 1.e4)
    flux_r = flux_r.clip(0.001, 1.e4)
    flux_i = flux_i.clip(0.001, 1.e4)

    # clip ivar as well, although this should not really be a problem as ivar
    # seems to always be well behaved

    ivar_g = ivar_g.clip(0.0, 70)
    ivar_r = ivar_r.clip(0.0, 70)
    ivar_i = ivar_i.clip(0.0, 70)

    ivar = ivar_g + ivar_r + ivar_i
    ivarsum = ivar_g + ivar_r + ivar_i

    fsum = flux_g*ivar_g + flux_r*ivar_r + flux_i*ivar_i

    flux = fsum/ivarsum

    return flux, ivarsum


def calc_c(modelflux, psfflux):
    return 1.0-psfflux/modelflux
