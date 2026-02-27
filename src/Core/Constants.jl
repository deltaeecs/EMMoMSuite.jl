module Constants

export c0, mu0, eps0, eta0

"Speed of light in vacuum (m/s)"
const c0 = 299792458.0

"Permeability of free space (H/m)"
const mu0 = 4π * 1e-7

"Permittivity of free space (F/m)"
const eps0 = 1.0 / (c0^2 * mu0)

"Intrinsic impedance of free space (Ohm)"
const eta0 = sqrt(mu0 / eps0)

end
