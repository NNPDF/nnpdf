#!awk -f
BEGIN { s = 38.8*38.8 } 
NR > 1 { 
	d = $8				# Data
	E = sqrt($9*$9+$10*$10)		# Error
	h = $4-$2                       # Mass bin-width

	# Bin-center values
	xFc = $1			  # Feynman-x
	Mc = h/2.0 + $2			  # Invariant mass
	Jc = sqrt(xFc*xFc + 4.0*Mc*Mc/s)  # Jacobian

	# Bin-average values
	xFa = $6
	Ma = $5
	Ja = sqrt(xFa*xFa + 4.0*Ma*Ma/s)

	yc  = 0.5*log(( Jc + xFc )/( Jc - xFc))
	ya  = 0.5*log(( Ja + xFa )/( Ja - xFa))

	print Ja, 0 
}
END{}     # End section
