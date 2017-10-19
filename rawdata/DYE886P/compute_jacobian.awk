#!awk -f
BEGIN { s = 38.8*38.8 } 
NR > 1 { 
	d = $7				# Data
	E = sqrt($8*$8+$9*$9)		# Error
	h = $3-$2                       # Mass bin-width

	# Bin-center values
	xFc = $1			  # Feynman-x
	Mc = h/2.0 + $2			  # Invariant mass
	Jc = sqrt(xFc*xFc + 4.0*Mc*Mc/s)  # Jacobian

	# Bin-average values
	xFa = $5
	Ma = $4
	Ja = sqrt(xFa*xFa + 4.0*Ma*Ma/s)

	yc  = 0.5*log(( Jc + xFc )/( Jc - xFc))
	ya  = 0.5*log(( Ja + xFa )/( Ja - xFa))

	print Ja
}
END{}     # End section
