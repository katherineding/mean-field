#Some mean field notes from discussion with Tom 09/09/2019:

In general, Mean field Hamiltonian does not respect all the symmetries the canonical hubbard hamiltonian has (such as particle number, spin conservation) because we did partial contractions and changed the physics of system 

We choose a basis for the MF hamiltonian that is consistent with its symmetries (this changes depending on how we define the mean fields). This basis would not be a complete basis for the hubbard hilbert space, but is instead only dimension O(N), single "particle"/'quasiparticle' hilbert space, making it easy to handle. 

If we choose a (n_up, n_dn) basis, which is what we'd usually expect, with all daggers in one vector and all annihilation operators in another vector. This basis, or Hilbert subspace, does NOT contain superconductivity

If we choose a (particle, hole) basis, then we'd have $<ci_up ci_dn>$ terms. These brackets evaluate to numbers in this basis. This basis or Hilbert subspace CAN contain superconductivity

c_i_up^dagger <-------time reversal -----> c_i_dn

Code shows Stripe order competes with SC order. Stripes exist ==> no SC


Code result strongly depends on initialization, and convergence is not gauranteed. There are ways to "force" things to converge, like asjusting the chemical potential each iteration to make sure filling level stays correct, or instead of using new values completely, use a linear combination of old and new values. But those techniques are kinda iffy too

With Delta = 0, the up electrons and down electrons don't talk to each other except between iterations, when we update n and sigma and hence change the mean field interaction.

All in all, this is mean field, its not exactly correct no matter what tricks we use to make it behave better
===
Mathematica alternatives = Octave, Maple, Python