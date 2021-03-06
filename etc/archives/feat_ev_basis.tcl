# EV title
set fmri(evtitle$ev_num) "$ev_name"

# Basic waveform shape (EV $ev_num)
# 0 : Square
# 1 : Sinusoid
# 2 : Custom (1 entry per volume)
# 3 : Custom (3 column format)
# 4 : Interaction
# 10 : Empty (all zeros)
set fmri(shape$ev_num) 3

# Convolution (EV $ev_num)
# 0 : None
# 1 : Gaussian
# 2 : Gamma
# 3 : Double-Gamma HRF
# 4 : Gamma basis functions
# 5 : Sine basis functions
# 6 : FIR basis functions
set fmri(convolve$ev_num) 7

# Convolve phase (EV $ev_num)
set fmri(convolve_phase$ev_num) 0

# Apply temporal filtering (EV $ev_num)
set fmri(tempfilt_yn$ev_num) $tempfilt_yn

# Add temporal derivative (EV $ev_num)
set fmri(deriv_yn$ev_num) $temporalderiv

# Custom EV file (EV $ev_num)
set fmri(custom$ev_num) "$cond_file"

# FIR basis functions number (EV 1)
set fmri(basisfnum1) 3

# Optimal/custom HRF convolution file (EV 1)
set fmri(bfcustom1) "$basis_file"	# $FSLDIR/etc/default_flobs.flobs/hrfbasisfns.txt

# Orth basis functions wrt each other
set fmri(basisorth1) $basis_orth	# 0 or 1
