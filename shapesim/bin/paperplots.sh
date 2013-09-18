for type in T flux; do

    if [[ type == "T" ]]; then
        panel='(a)'
    else
        panel='(b)'
    fi
    python plot-cshapesim.py \
        --s2n-field "${type}_s2n" \
        cbafit-geg02rcomb,cbafit-geg06rcomb,cbafit-deg01rcomb,cbafit-deg02rcomb 0.01 \
        -y -0.021,0.021 \
        --labels "exp \sigma_g/\sigma_{PSF}=2.0","exp \sigma_g/\sigma_{PSF}=1.4","dev \sigma_g/\sigma_{PSF}=2.0","dev \sigma_g/\sigma_{PSF}=1.4" \
        --eps "/astro/u/esheldon/tmp/cbafit-geg-${type}-s2n.eps" \
        --panel "$panel"
done

