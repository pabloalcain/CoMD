# CoMD

Implementation of the Constrained Molecular Dynamics, based on
[Bonasera et al. PRC64 024612 (2001)][1]. Some changes are: we won't
add Pauli blocking, and not necessarily check in every step whether f
gets larger than 1.

The potentials and the distribution functions are all integrated, so
the "gaussian" characteristic of the wavepacket is embedded on the
expression for the interaction and the occupation probability.


[1]: http://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.024612

## Usage

To compile, go to src/ and run `make`. It will create an executable
called comd.x

So far this executable comes from a *barebone* main function where you
can choose how many timesteps to run on runtime, with `./comd.x
<nsteps>` but everything else has to be tweaked from source. It will
be changed either to get everything from input files or exposing
classes to python through boost.