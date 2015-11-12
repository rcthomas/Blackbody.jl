
using Planck
using Base.Test

@test_approx_eq planck( [ 5000.0 ], 5000.0, PerHertz() ) [ 1.00964e-5 ]
