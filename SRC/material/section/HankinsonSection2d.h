/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** ****************************************************************** */

// Description: This file contains the class definition for HankinsonSection2d.
//
// HankinsonSection2d is a 2-DOF coupled SectionForceDeformation for the
// M-N-V interaction of timber dowel-type connections, following the
// modelling approach of Schweigler et al.: a dowel's in-plane slip
// (d_parallel, d_perpendicular) is fed to a scalar slip law that is
// interpolated between a 0-degree and a 90-degree calibration by the
// Hankinson criterion (exponent n = 2),
//
//     X_alpha = X_0 * X_90 / (X_0 * sin^2(alpha) + X_90 * cos^2(alpha))
//
// where alpha is the (evolving) load-to-grain angle of the current slip
// vector. This is the literal Schweigler formulation: a single
// UniaxialMaterial's named parameters (via the standard setParameter/
// updateParameter mechanism) are Hankinson-interpolated between a 0-deg and
// a 90-deg value and pushed into the material at every commit -- see
// DowelType and SAWSMaterial, which were both extended with setParameter/
// updateParameter for exactly this purpose.
//
// An earlier version of this class also offered a response-level blending
// mode (two independently-calibrated UniaxialMaterial copies, driven by the
// same signed scalar slip, with their scalar force responses Hankinson-
// blended). That mode was removed: blending the *responses* of two
// materials that each carry their own independent hysteretic history has no
// clear physical basis -- the Hankinson criterion is a statement about how
// a single physical law's own parameters vary with grain angle, not about
// how to combine two unrelated internal states. The parameter-level mode
// below is the only one implemented; it requires the driven material to
// expose the relevant parameters via setParameter/updateParameter, which is
// a small, self-contained addition for any given UniaxialMaterial (see
// DowelType::setParameter / SAWSMaterial::setParameter for two examples).
//
// This section is meant to be used with the existing ZeroLengthSection
// element (order-2, codes {SECTION_RESPONSE_P, SECTION_RESPONSE_VY} by
// default), which builds a fully consistent coupled element stiffness
// K = A^T * kb * A from whatever 2x2 tangent this section returns -- no new
// element is required.
//
// See connection-mn-interaction/schweigler_implementation for the
// reference (monotonic-only) Python implementation this generalizes.

#ifndef HankinsonSection2d_h
#define HankinsonSection2d_h

#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Parameter.h>

#include <vector>
#include <string>

class UniaxialMaterial;

class HankinsonSection2d : public SectionForceDeformation
{
  public:
    enum SignRule  { SIGN_ROTATING = 0, SIGN_FIXED = 1 };

    // Parameter-level Hankinson interpolation -- one material whose named
    // parameters are re-interpolated between 0-deg and 90-deg values and
    // pushed into it (via setParameter/updateParameter) at every commit.
    HankinsonSection2d(int tag,
                      UniaxialMaterial &mat,
                      int nPar, const char **parNames,
                      const double *val0, const double *val90,
                      int code1 = SECTION_RESPONSE_P,
                      int code2 = SECTION_RESPONSE_VY,
                      int signRule = SIGN_ROTATING,
                      bool resync = true,
                      double gcapFac = 10.0,
                      double dtolRel = 1.0e-10);

    HankinsonSection2d();          // blank object for FEM_ObjectBroker
    ~HankinsonSection2d();

    const char *getClassType(void) const { return "HankinsonSection2d"; }

    int setTrialSectionDeformation(const Vector &e);
    const Vector &getSectionDeformation(void) { return e_t; }
    const Vector &getStressResultant(void)    { return q_t; }
    const Matrix &getSectionTangent(void)     { return kb_t; }
    const Matrix &getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    SectionForceDeformation *getCopy(void);
    const ID &getType(void);
    int getOrder(void) const { return 2; }

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &info);
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:

  private:
    // ---- configuration (fixed for the life of the object) -------------
    int signRule;       // SignRule
    int code1, code2;
    double gcapFac;      // tangent radial-stiffness cap, relative to max(E_par,E_perp)
    double dtolRel;      // slip-direction tolerance, relative to CdeltaMax
    bool resyncAfterMutation;

    UniaxialMaterial *theMat;       // the single driven material
    double E_par, E_perp;           // cached initial tangents (0-deg / 90-deg)

    // ---- parameter table -----------------------------------
    int numPar;
    std::vector<std::string> parNames;
    std::vector<double> parVal0, parVal90;   // as given, for getCopy()/Print()/sendSelf()
    std::vector<double> parA, parB;          // 1/X0 , 1/X90 - 1/X0
    std::vector<Parameter> thePar;
    void setupParameters(void);       // binds parNames[i] -> thePar[i] on theMat
    void updateLaws(double w);        // mutate + (optionally) resync; commit/revert ONLY

    // ---- trial state ------------------------------------------------
    Vector e_t, q_t;
    Matrix kb_t;
    double n0_t, n1_t;     // unit slip direction actually used this step
    double delta_t;        // ||e||
    double s_t;             // signed scalar slip fed to the sub-material
    double w_t;              // sin^2(alpha) used this step
    double F_t, k_t;          // scalar force / tangent from evaluateLaw
    int    sigma_t;
    bool   virgin_t;

    // ---- committed state ----------------------------------------------
    Vector Ce, Cq;
    Matrix Ckb;
    double Cn0, Cn1;        // committed reference slip axis (0,0) == virgin
    int    Csigma;
    double Cdelta, Cs, Cw, CF;
    double CdeltaMax;       // running max committed delta, for the direction tolerance
    double CexplicitErr;    // staggering-error monitor (relative force jump at commit)

    Matrix K0Mat;            // scratch for getInitialTangent()
    ID     theCode;          // scratch for getType()

    void evaluateLaw(double s, double &F, double &k);
    void buildInitialTangent(void);   // computes E_par / E_perp
};

#endif
