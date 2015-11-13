using Base.Test

using Blackbody

for i = 1:20
    temperature = 1000.0 * i
    @test_approx_eq_eps quadgk( wl -> planck( wl, temperature, PerAngstrom ), 0.0, 1.0e9 )[ 1 ] / stefan_boltzmann( temperature ) 1.0 1.0e-7
end
