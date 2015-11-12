module Blackbody

using PhysicalConstants

export PerHertz, PerAngstrom, planck, wien

const c0 = PhysicalConstants.CGS.PlancksConstantH * PhysicalConstants.CGS.SpeedOfLight
const c1 = 2.0 * c0
const c2 = c0 / PhysicalConstants.CGS.Boltzmann
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
    [ planck( 1.0e8 / wavelength[ i ], temperature, output ) for i = 1:length( wavelength ) ]
end

function planck( x::Number, temperature::Number, output::Type{PerHertz} )
    c1 * x ^ 3 / expm( c2 * x / temperature )
end

function planck( x::Number, temperature::Number, output::Type{PerAngstrom} )
    c1 * x ^ 5 / expm( c2 * x / temperature )
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