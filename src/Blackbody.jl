module Blackbody

using PhysicalConstants

export PerHertz, PerAngstrom, planck, wien

const c0  = PhysicalConstants.CGS.PlancksConstantH * PhysicalConstants.CGS.SpeedOfLight
const c1h = 2.0 * c0 * PhysicalConstants.CGS.SpeedOfLight
const c1a = c1h / 1.0e8
const c2  = c0 / PhysicalConstants.CGS.Boltzmann

const b1 = 1.0e8 * c2 / 2.821439372
const b2 = 1.0e8 * c2 / 4.965114231744276

abstract PerUnit

immutable PerHertz <: PerUnit
end

immutable PerAngstrom <: PerUnit
end

"""
Planck (blackbody) spectral radiance or specific intensity in CGS units.

    Args:
        wavelength (Array{Number,1} or Number): Wavelength(s) in Angstroms.
        temperature (Number): Temperature in Kelvin.
        output (Type{PerUnit}): Output convention (PerHertz or PerAngstrom).
    
    Returns:
        Array{Float64,1} or Float64: Planck spectral radiance in CGS units.
"""
function planck{T1<:Number,T2<:PerUnit}( wavelength::Array{T1,1}, temperature::Number, output::Type{T2} )
    [ planck( wavelength[ i ], temperature, output ) for i = 1:length( wavelength ) ]
end

function planck( wavelength::Number, temperature::Number, output::Type{PerHertz} )
    x = 1.0e8 / wavelength
    c1h * x ^ 3 / expm1( c2 * x / temperature )
end

function planck( wavelength::Number, temperature::Number, output::Type{PerAngstrom} )
    x = 1.0e8 / wavelength
    c1a * x ^ 5 / expm1( c2 * x / temperature )
end

"""
Wavelength at the maximum Planck spectral radiance for a given temperature
(Wien's displacement law).

    Args:
        temperature (Number): Temperature in Kelvin.
        output (Type{PerUnit}): Unit convention (PerHertz or PerAngstrom).
    
    Returns:
        Float64: Wavelength at maximum Planck spectral radiance in Angstroms.
"""
wien( temperature::Number, output::Type{PerHertz} ) = b1 / temperature

wien( temperature::Number, output::Type{PerAngstrom} ) = b2 / temperature

end

# using Blackbody
#
# wavelengths = collect( logspace( 3, 4, 100000000 ) )
# for ( wl, bb ) in zip( wavelengths, planck( wavelengths, 10000.0, PerAngstrom ) )
#     println( "$wl $bb" )
# end
# println()
# 
# wlmax = wien( 10000.0, PerAngstrom )
# println( wlmax, " ", planck( wlmax, 10000.0, PerAngstrom ) )
# println() 
# 
# @time planck( wavelengths, 10000.0, PerAngstrom )
# @time planck( wavelengths, 10000.0, PerAngstrom )
# @time planck( wavelengths, 10000.0, PerAngstrom )
# @time planck( wavelengths, 10000.0, PerAngstrom )
# @time planck( wavelengths, 10000.0, PerAngstrom )


