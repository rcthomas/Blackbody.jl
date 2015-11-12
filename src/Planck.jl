
module Planck

using PhysicalConstants

export PerHertz, PerAngstrom, planck

const c0 = PhysicalConstants.CGS.PlancksConstantH * PhysicalConstants.CGS.SpeedOfLight
const c1 = 2.0 * c0
const c2 = c0 / PhysicalConstants.CGS.Boltzmann
const b  = 1.0e8 * c2 / 4.965114231744276

abstract PerUnit

immutable PerHertz <: PerUnit
end

immutable PerAngstrom <: PerUnit
end

"""
Planck (blackbody) spectral radiance or specific intensity in CGS units.

    Args:
        wavelength (Array{Number,1}): Wavelengths in Angstroms.
        temperature (Number): Temperature in Kelvin.
        output (PerUnit): Output convention (PerHertz or PerAngstrom).
    
    Returns:
        Array{Float64,1}: Spectral radiance in erg cm^-2 sr^-1 Hz^-1 or AA^-1.
"""
function planck{T<:Number}( wavelength::Array{T,1}, temperature::Number, output::PerUnit )
    impl( 1.0e8 ./ wavelength, temperature, output )
end

function impl{T<:Number}( x::Array{T,1}, temperature::Number, output::PerHertz )
    c1 .* x .^ 3 ./ expm1( c2 .* x ./ temperature )
end

function impl{T<:Number}( x::Array{T,1}, temperature::Number, output::PerAngstrom )
    c1 .* x .^ 5 ./ expm1( c2 .* x ./ temperature )
end

end
