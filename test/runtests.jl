
using Blackbody
using Base.Test

#@test_approx_eq planck( [ 5000.0 ], 5000.0, PerHertz() ) [ 1.00964e-5 ] make better test

f( wl ) = planck( wl, 12000.0, PerAngstrom )

@test_approx_eq_eps quadgk( f, 0.0, 1.0e9 )[ 1 ] / stefan_boltzmann( 12000.0 ) 1.0 1.0e-14


