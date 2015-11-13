using Base.Test

using Blackbody

custom_handler( r::Test.Success ) = println( "PASS $(r.expr)" )
custom_handler( r::Test.Failure ) = println( "FAIL $(r.expr)" )
custom_handler( r::Test.Error   ) = rethrow( r )

Test.with_handler( custom_handler ) do
    for i = 1:20
        temperature = 1000.0 * i
        @test_approx_eq_eps quadgk( wl -> planck( wl, temperature, PerAngstrom ), 0.0, 1.0e9 )[ 1 ] / stefan_boltzmann( temperature ) 1.0 1.0e-7
    end
end
