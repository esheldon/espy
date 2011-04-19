import es_sdsspy
s82=es_sdsspy.sg.Stripe82Epochs()

s82.collate()
s82.extract_psf_fwhm()
s82.collate_psf_fwhm()
