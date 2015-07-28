#ifndef POTENTIAL_H
#define POTENTIAL_H

class Potential
{
    /* Basic CoMD potential structure, with all the parameters */
    double t0, r0;
    double t3;
    double asym;
    double cs;
    double e;
    
public:
    double rcut;
    Potential();
    double dphi(double dr);
};

#endif
