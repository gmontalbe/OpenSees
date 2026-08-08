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

// See HankinsonSection2d.h for a full description.
//
// Tcl / OpenSeesPy command:
//
//   section HankinsonSection2d $secTag -paramHankinson $matTag \
//       -param F0 19.1 14.18  -param K0 24460. 19640.  -param K1 462. 739. \
//       <-codes P Vy> <-signRule rotating|fixed> <-resync 1> <-gcap 10.0> <-dtol 1e-10>

#include <HankinsonSection2d.h>

#include <UniaxialMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Response.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

#include <math.h>
#include <float.h>
#include <string.h>

// ---------------------------------------------------------------------
// Tcl/Python parser
// ---------------------------------------------------------------------

static int
HD2D_parseCode(const char *s)
{
    if (strcmp(s, "P") == 0)  return SECTION_RESPONSE_P;
    if (strcmp(s, "Vy") == 0) return SECTION_RESPONSE_VY;
    if (strcmp(s, "My") == 0) return SECTION_RESPONSE_MY;
    if (strcmp(s, "Mz") == 0) return SECTION_RESPONSE_MZ;
    if (strcmp(s, "Vz") == 0) return SECTION_RESPONSE_VZ;
    if (strcmp(s, "T") == 0)  return SECTION_RESPONSE_T;
    return SECTION_RESPONSE_NONE;
}

void *
OPS_HankinsonSection2d()
{
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 4) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: section HankinsonSection2d tag -paramHankinson mat? -param name v0 v90 ...\n";
        return 0;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
        opserr << "WARNING invalid section HankinsonSection2d tag" << endln;
        return 0;
    }

    const char *modeFlag = OPS_GetString();
    if (strcmp(modeFlag, "-paramHankinson") != 0) {
        opserr << "WARNING unknown option \"" << modeFlag
               << "\" for section HankinsonSection2d " << tag
               << " -- expected -paramHankinson" << endln;
        return 0;
    }

    int matTag;
    numData = 1;
    if (OPS_GetIntInput(&numData, &matTag) != 0) {
        opserr << "WARNING invalid mat tag for section HankinsonSection2d " << tag << endln;
        return 0;
    }
    UniaxialMaterial *mat = OPS_getUniaxialMaterial(matTag);
    if (mat == 0) {
        opserr << "WARNING uniaxialMaterial " << matTag
               << " not found for section HankinsonSection2d " << tag << endln;
        return 0;
    }

    // ---- optional flags -------------------------------------------------
    int code1 = SECTION_RESPONSE_P, code2 = SECTION_RESPONSE_VY;
    int signRule = HankinsonSection2d::SIGN_ROTATING;
    double gcap = 10.0, dtol = 1.0e-10;
    bool resync = true;

    std::vector<std::string> parNames;
    std::vector<double> val0, val90;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *flag = OPS_GetString();

        if (strcmp(flag, "-codes") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 2) {
                opserr << "WARNING -codes needs two response codes" << endln;
                return 0;
            }
            const char *s1 = OPS_GetString();
            const char *s2 = OPS_GetString();
            code1 = HD2D_parseCode(s1);
            code2 = HD2D_parseCode(s2);
            if (code1 == SECTION_RESPONSE_NONE || code2 == SECTION_RESPONSE_NONE) {
                opserr << "WARNING invalid -codes \"" << s1 << "\" \"" << s2 << "\"" << endln;
                return 0;
            }
        } else if (strcmp(flag, "-signRule") == 0) {
            const char *s = OPS_GetString();
            if (strcmp(s, "rotating") == 0)
                signRule = HankinsonSection2d::SIGN_ROTATING;
            else if (strcmp(s, "fixed") == 0)
                signRule = HankinsonSection2d::SIGN_FIXED;
            else {
                opserr << "WARNING invalid -signRule \"" << s << "\" (want rotating|fixed)" << endln;
                return 0;
            }
        } else if (strcmp(flag, "-gcap") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &gcap) != 0) {
                opserr << "WARNING invalid -gcap value" << endln;
                return 0;
            }
        } else if (strcmp(flag, "-dtol") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &dtol) != 0) {
                opserr << "WARNING invalid -dtol value" << endln;
                return 0;
            }
        } else if (strcmp(flag, "-resync") == 0) {
            int iv;
            numData = 1;
            if (OPS_GetIntInput(&numData, &iv) != 0) {
                opserr << "WARNING invalid -resync value" << endln;
                return 0;
            }
            resync = (iv != 0);
        } else if (strcmp(flag, "-param") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING -param needs name, val0 and val90" << endln;
                return 0;
            }
            const char *name = OPS_GetString();
            double v[2];
            numData = 2;
            if (OPS_GetDoubleInput(&numData, v) != 0) {
                opserr << "WARNING invalid -param " << name << " values" << endln;
                return 0;
            }
            // Hankinson's formula only needs val0 and val90 to be nonzero and
            // same-signed (e.g. DowelType's "negative-side" parameters such as
            // F0n/Dcn are, by convention, negative magnitudes -- the formula
            // is sign-invariant: interpolating two same-signed values always
            // returns a value of that sign).
            if (v[0] * v[1] <= 0.0) {
                opserr << "WARNING -param " << name
                       << ": Hankinson interpolation requires val0 and val90 to be "
                       << "nonzero and of the same sign, got "
                       << v[0] << " and " << v[1] << endln;
                return 0;
            }
            parNames.push_back(std::string(name));
            val0.push_back(v[0]);
            val90.push_back(v[1]);
        } else {
            opserr << "WARNING unknown option \"" << flag
                   << "\" for section HankinsonSection2d " << tag << endln;
            return 0;
        }
    }

    if (parNames.empty()) {
        opserr << "WARNING section HankinsonSection2d " << tag
               << " -paramHankinson requires at least one -param entry" << endln;
        return 0;
    }
    std::vector<const char *> nameptrs(parNames.size());
    for (size_t i = 0; i < parNames.size(); i++)
        nameptrs[i] = parNames[i].c_str();

    return new HankinsonSection2d(tag, *mat, (int)parNames.size(), &nameptrs[0],
                                 &val0[0], &val90[0], code1, code2,
                                 signRule, resync, gcap, dtol);
}

// ---------------------------------------------------------------------
// Constructors / destructor
// ---------------------------------------------------------------------

HankinsonSection2d::HankinsonSection2d(int tag, UniaxialMaterial &mat,
                                    int nPar, const char **names,
                                    const double *val0, const double *val90,
                                    int c1, int c2, int sRule, bool resync,
                                    double gcapF, double dtolR)
: SectionForceDeformation(tag, SEC_TAG_HankinsonSection2d),
  signRule(sRule),
  code1(c1), code2(c2), gcapFac(gcapF), dtolRel(dtolR),
  resyncAfterMutation(resync),
  theMat(0),
  E_par(0.0), E_perp(0.0),
  numPar(nPar),
  parNames(nPar), parVal0(nPar), parVal90(nPar), parA(nPar), parB(nPar), thePar(nPar),
  e_t(2), q_t(2), kb_t(2, 2),
  n0_t(0.0), n1_t(0.0), delta_t(0.0), s_t(0.0), w_t(0.0), F_t(0.0), k_t(0.0),
  sigma_t(1), virgin_t(true),
  Ce(2), Cq(2), Ckb(2, 2),
  Cn0(0.0), Cn1(0.0), Csigma(1),
  Cdelta(0.0), Cs(0.0), Cw(0.0), CF(0.0),
  CdeltaMax(0.0), CexplicitErr(0.0),
  K0Mat(2, 2), theCode(2)
{
    theMat = mat.getCopy();
    if (theMat == 0) {
        opserr << "HankinsonSection2d::HankinsonSection2d -- failed to get a copy of the material\n";
        exit(-1);
    }

    for (int i = 0; i < numPar; i++) {
        parNames[i] = names[i];
        parVal0[i]  = val0[i];
        parVal90[i] = val90[i];
        parA[i] = 1.0 / val0[i];
        parB[i] = 1.0 / val90[i] - 1.0 / val0[i];
    }
    this->setupParameters();

    // Prime the sub-material at both calibration angles (harmless -- no
    // history has accumulated yet) to get E_par/E_perp, then leave it at
    // the 0-deg (alpha = 0) parameterization.
    this->updateLaws(0.0);
    E_par = theMat->getInitialTangent();
    this->updateLaws(1.0);
    E_perp = theMat->getInitialTangent();
    this->updateLaws(0.0);
    theMat->revertToStart();

    this->revertToStart();
}

HankinsonSection2d::HankinsonSection2d()
: SectionForceDeformation(0, SEC_TAG_HankinsonSection2d),
  signRule(SIGN_ROTATING),
  code1(SECTION_RESPONSE_P), code2(SECTION_RESPONSE_VY),
  gcapFac(10.0), dtolRel(1.0e-10),
  resyncAfterMutation(true),
  theMat(0),
  E_par(0.0), E_perp(0.0),
  numPar(0),
  e_t(2), q_t(2), kb_t(2, 2),
  n0_t(0.0), n1_t(0.0), delta_t(0.0), s_t(0.0), w_t(0.0), F_t(0.0), k_t(0.0),
  sigma_t(1), virgin_t(true),
  Ce(2), Cq(2), Ckb(2, 2),
  Cn0(0.0), Cn1(0.0), Csigma(1),
  Cdelta(0.0), Cs(0.0), Cw(0.0), CF(0.0),
  CdeltaMax(0.0), CexplicitErr(0.0),
  K0Mat(2, 2), theCode(2)
{
    theMat = 0;
}

HankinsonSection2d::~HankinsonSection2d()
{
    if (theMat != 0) delete theMat;
}

void
HankinsonSection2d::buildInitialTangent(void)
{
    // E_par/E_perp are computed in the constructor, before setupParameters()/
    // updateLaws() are meaningful to call from here; nothing to do here.
}

// ---------------------------------------------------------------------
// Parameter binding and mutation
// ---------------------------------------------------------------------

void
HankinsonSection2d::setupParameters(void)
{
    for (int i = 0; i < numPar; i++) {
        thePar[i] = Parameter();
        const char *argv[1];
        argv[0] = parNames[i].c_str();
        int ok = theMat->setParameter(argv, 1, thePar[i]);
        if (ok < 0) {
            opserr << "HankinsonSection2d::setupParameters -- material \""
                   << theMat->getClassType() << "\" does not support parameter \""
                   << parNames[i].c_str() << "\"" << endln;
        }
    }
}

// Recompute X_alpha for every bound parameter and push it into the
// sub-material. This mutates the sub-material's envelope, so it must ONLY
// be called from commitState()/revertToLastCommit()/revertToStart() -- see
// HankinsonSection2d.h and the ASDCoupledHinge3D precedent this follows.
void
HankinsonSection2d::updateLaws(double w)
{
    for (int i = 0; i < numPar; i++) {
        double denom = parA[i] + parB[i] * w;
        double Xw = (denom != 0.0) ? 1.0 / denom : parVal0[i];
        thePar[i].update(Xw);
    }

    if (resyncAfterMutation) {
        // Generic "envelope changed underneath you" re-sync recipe (see
        // Pinching4Material::updateParameter): recompute the stress at the
        // strain we are already at, and commit that as the current state.
        theMat->setTrialStrain(theMat->getStrain());
        theMat->commitState();
    }
}

// ---------------------------------------------------------------------
// Core kinematics / scalar law evaluation
// ---------------------------------------------------------------------

void
HankinsonSection2d::evaluateLaw(double s, double &F, double &k)
{
    theMat->setTrialStrain(s);
    F = theMat->getStress();
    k = theMat->getTangent();
}

int
HankinsonSection2d::setTrialSectionDeformation(const Vector &e)
{
    e_t(0) = e(0);
    e_t(1) = e(1);
    double e0 = e(0), e1 = e(1);
    delta_t = sqrt(e0 * e0 + e1 * e1);

    double dtol = dtolRel * CdeltaMax;
    double dtolFloor = 1.0e-12 * fmax(1.0, CdeltaMax);
    if (dtol < dtolFloor) dtol = dtolFloor;

    if (delta_t > dtol) {
        n0_t = e0 / delta_t;
        n1_t = e1 / delta_t;
        if (Cn0 == 0.0 && Cn1 == 0.0)
            sigma_t = 1;
        else if (signRule == SIGN_FIXED)
            sigma_t = Csigma;
        else
            sigma_t = (e0 * Cn0 + e1 * Cn1 >= 0.0) ? 1 : -1;
        virgin_t = false;
    } else if (Cn0 != 0.0 || Cn1 != 0.0) {
        // Slip has returned to (near) zero on a previously-loaded axis:
        // hold the last committed direction so the residual force stays
        // collinear with it instead of collapsing to an undefined direction.
        // Cn0/Cn1 IS the stable "positive-s" direction n-hat (= sigma*n from
        // the last real measurement, see commitState()), so sigma_t must be
        // +1 here -- NOT Csigma, which would double-apply the sign (Cn0 is
        // already sigma-adjusted) and silently flip the reference axis on
        // the next commit. That flip is invisible for an odd force law
        // (F(-s)=-F(s), e.g. any purely elastic material) because the two
        // sign errors cancel in q = sigma*F(sigma*delta)*n, but corrupts the
        // sign for any material with a residual/asymmetric force.
        n0_t = Cn0;
        n1_t = Cn1;
        sigma_t = 1;
        virgin_t = false;
    } else {
        // True virgin state: no direction has ever been established.
        virgin_t = true;
        n0_t = 0.0; n1_t = 0.0; sigma_t = 1;
        s_t = 0.0; w_t = Cw; F_t = 0.0; k_t = 0.5 * (E_par + E_perp);
        q_t.Zero();
        kb_t = this->getInitialTangent();
        theMat->setTrialStrain(0.0);
        return 0;
    }

    s_t = sigma_t * delta_t;
    w_t = Cw;   // committed (lagged) alpha -- see updateLaws()/commitState()

    this->evaluateLaw(s_t, F_t, k_t);

    q_t(0) = sigma_t * F_t * n0_t;
    q_t(1) = sigma_t * F_t * n1_t;

    double g = (s_t != 0.0) ? F_t / s_t : 0.5 * (E_par + E_perp);
    double gcap = gcapFac * fmax(E_par, E_perp);
    if (g >  gcap) g =  gcap;
    if (g < -gcap) g = -gcap;

    double nn00 = n0_t * n0_t, nn01 = n0_t * n1_t, nn11 = n1_t * n1_t;
    kb_t(0, 0) = k_t * nn00 + g * (1.0 - nn00);
    kb_t(0, 1) = (k_t - g) * nn01;
    kb_t(1, 0) = kb_t(0, 1);
    kb_t(1, 1) = k_t * nn11 + g * (1.0 - nn11);

    return 0;
}

const Matrix &
HankinsonSection2d::getInitialTangent(void)
{
    K0Mat(0, 0) = E_par;
    K0Mat(0, 1) = 0.0;
    K0Mat(1, 0) = 0.0;
    K0Mat(1, 1) = E_perp;
    return K0Mat;
}

// ---------------------------------------------------------------------
// State management
// ---------------------------------------------------------------------

int
HankinsonSection2d::commitState(void)
{
    int err = theMat->commitState();
    if (err != 0)
        return err;

    if (!virgin_t) {
        Cn0 = sigma_t * n0_t;
        Cn1 = sigma_t * n1_t;
        Csigma = sigma_t;
    }
    Ce = e_t; Cq = q_t; Ckb = kb_t;
    Cdelta = delta_t; Cs = s_t; CF = F_t;
    if (delta_t > CdeltaMax) CdeltaMax = delta_t;

    double wNew = n1_t * n1_t;
    double Fbefore = theMat->getStress();
    this->updateLaws(wNew);
    double Fafter = theMat->getStress();
    CexplicitErr = fabs(Fafter - Fbefore) / fmax(fabs(Fbefore), 1.0);
    Cw = wNew;

    return 0;
}

int
HankinsonSection2d::revertToLastCommit(void)
{
    int err = theMat->revertToLastCommit();

    e_t = Ce; q_t = Cq; kb_t = Ckb;
    delta_t = Cdelta; s_t = Cs; F_t = CF; w_t = Cw;
    sigma_t = Csigma; n0_t = Cn0; n1_t = Cn1;
    virgin_t = (Cn0 == 0.0 && Cn1 == 0.0 && Cdelta <= 0.0);

    this->updateLaws(Cw);

    return err;
}

int
HankinsonSection2d::revertToStart(void)
{
    int err = theMat->revertToStart();

    Ce.Zero(); Cq.Zero();
    Ckb = this->getInitialTangent();
    Cn0 = 0.0; Cn1 = 0.0; Csigma = 1;
    Cdelta = 0.0; Cs = 0.0; Cw = 0.0; CF = 0.0;
    CdeltaMax = 0.0; CexplicitErr = 0.0;

    this->updateLaws(0.0);

    return err + this->revertToLastCommit();
}

// ---------------------------------------------------------------------
// Copy / type / serialization
// ---------------------------------------------------------------------

SectionForceDeformation *
HankinsonSection2d::getCopy(void)
{
    std::vector<const char *> nameptrs(numPar);
    for (int i = 0; i < numPar; i++)
        nameptrs[i] = parNames[i].c_str();

    // getCopy() must route through a real constructor (not a memberwise
    // copy): setupParameters() re-binds fresh Parameters to a fresh
    // getCopy() of theMat. A memberwise copy would leave every element's
    // thePar[i] pointing at the PROTOTYPE's material, so all elements
    // sharing this section tag would silently mutate one shared law.
    HankinsonSection2d *theCopy =
        new HankinsonSection2d(this->getTag(), *theMat,
                              numPar, numPar > 0 ? &nameptrs[0] : 0,
                              numPar > 0 ? &parVal0[0] : 0,
                              numPar > 0 ? &parVal90[0] : 0,
                              code1, code2, signRule, resyncAfterMutation,
                              gcapFac, dtolRel);

    theCopy->Ce = Ce; theCopy->Cq = Cq; theCopy->Ckb = Ckb;
    theCopy->Cn0 = Cn0; theCopy->Cn1 = Cn1; theCopy->Csigma = Csigma;
    theCopy->Cdelta = Cdelta; theCopy->Cs = Cs; theCopy->Cw = Cw; theCopy->CF = CF;
    theCopy->CdeltaMax = CdeltaMax;
    theCopy->CexplicitErr = CexplicitErr;
    theCopy->revertToLastCommit();

    return theCopy;
}

const ID &
HankinsonSection2d::getType(void)
{
    theCode(0) = code1;
    theCode(1) = code2;
    return theCode;
}

int
HankinsonSection2d::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    int dbTag = this->getDbTag();

    ID idData(6);
    idData(0) = this->getTag();
    idData(1) = signRule;
    idData(2) = code1;
    idData(3) = code2;
    idData(4) = resyncAfterMutation ? 1 : 0;
    idData(5) = numPar;
    res = theChannel.sendID(dbTag, cTag, idData);
    if (res < 0) {
        opserr << "HankinsonSection2d::sendSelf -- failed to send ID data\n";
        return res;
    }

    ID matData(2);
    matData(0) = theMat->getClassTag();
    int matDbTag = theMat->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        if (matDbTag != 0)
            theMat->setDbTag(matDbTag);
    }
    matData(1) = matDbTag;
    res = theChannel.sendID(dbTag, cTag, matData);
    if (res < 0) {
        opserr << "HankinsonSection2d::sendSelf -- failed to send material ID data\n";
        return res;
    }

    Vector data(8 + 4 * numPar);
    data(0) = gcapFac; data(1) = dtolRel;
    data(2) = E_par; data(3) = E_perp;
    data(4) = Cn0; data(5) = Cn1; data(6) = Csigma;
    data(7) = Cdelta;
    for (int i = 0; i < numPar; i++) {
        data(8 + 4 * i)     = parVal0[i];
        data(8 + 4 * i + 1) = parVal90[i];
        data(8 + 4 * i + 2) = parA[i];
        data(8 + 4 * i + 3) = parB[i];
    }
    res = theChannel.sendVector(dbTag, cTag, data);
    if (res < 0) {
        opserr << "HankinsonSection2d::sendSelf -- failed to send Vector data\n";
        return res;
    }

    if (numPar > 0) {
        // pack the parameter names as one comma-separated string
        std::string joined;
        for (int i = 0; i < numPar; i++) {
            if (i > 0) joined += ",";
            joined += parNames[i];
        }
        ID lenData(1);
        lenData(0) = (int)joined.size();
        res = theChannel.sendID(dbTag, cTag, lenData);
        if (res < 0) return res;

        ID nameData((int)joined.size() > 0 ? (int)joined.size() : 1);
        for (size_t i = 0; i < joined.size(); i++)
            nameData((int)i) = (int)joined[i];
        res = theChannel.sendID(dbTag, cTag, nameData);
        if (res < 0) return res;
    }

    res = theMat->sendSelf(cTag, theChannel);
    if (res < 0) {
        opserr << "HankinsonSection2d::sendSelf -- failed to send material\n";
        return res;
    }

    return res;
}

int
HankinsonSection2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    int dbTag = this->getDbTag();

    ID idData(6);
    res = theChannel.recvID(dbTag, cTag, idData);
    if (res < 0) {
        opserr << "HankinsonSection2d::recvSelf -- failed to receive ID data\n";
        return res;
    }

    this->setTag(idData(0));
    signRule = idData(1);
    code1 = idData(2);
    code2 = idData(3);
    resyncAfterMutation = (idData(4) != 0);
    numPar = idData(5);

    ID matData(2);
    res = theChannel.recvID(dbTag, cTag, matData);
    if (res < 0) {
        opserr << "HankinsonSection2d::recvSelf -- failed to receive material ID data\n";
        return res;
    }
    int matClassTag = matData(0);
    int matDbTag = matData(1);

    Vector data(8 + 4 * numPar);
    res = theChannel.recvVector(dbTag, cTag, data);
    if (res < 0) {
        opserr << "HankinsonSection2d::recvSelf -- failed to receive Vector data\n";
        return res;
    }
    gcapFac = data(0); dtolRel = data(1);
    E_par = data(2); E_perp = data(3);
    Cn0 = data(4); Cn1 = data(5); Csigma = (int)data(6);
    Cdelta = data(7);

    parVal0.resize(numPar); parVal90.resize(numPar);
    parA.resize(numPar); parB.resize(numPar);
    parNames.resize(numPar); thePar.resize(numPar);
    for (int i = 0; i < numPar; i++) {
        parVal0[i] = data(8 + 4 * i);
        parVal90[i] = data(8 + 4 * i + 1);
        parA[i] = data(8 + 4 * i + 2);
        parB[i] = data(8 + 4 * i + 3);
    }

    if (numPar > 0) {
        ID lenData(1);
        res = theChannel.recvID(dbTag, cTag, lenData);
        if (res < 0) return res;
        int len = lenData(0);
        ID nameData(len > 0 ? len : 1);
        res = theChannel.recvID(dbTag, cTag, nameData);
        if (res < 0) return res;
        std::string joined;
        for (int i = 0; i < len; i++)
            joined += (char)nameData(i);
        // split on commas
        size_t start = 0;
        for (int i = 0; i < numPar; i++) {
            size_t comma = joined.find(',', start);
            if (comma == std::string::npos) comma = joined.size();
            parNames[i] = joined.substr(start, comma - start);
            start = comma + 1;
        }
    }

    if (theMat != 0) { delete theMat; theMat = 0; }

    theMat = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMat == 0) {
        opserr << "HankinsonSection2d::recvSelf -- could not get a UniaxialMaterial\n";
        return -1;
    }
    theMat->setDbTag(matDbTag);
    res = theMat->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
        opserr << "HankinsonSection2d::recvSelf -- material failed to recvSelf\n";
        return res;
    }

    this->setupParameters();
    this->revertToLastCommit();

    return 0;
}

void
HankinsonSection2d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "HankinsonSection2d, tag: " << this->getTag() << endln;
        s << "\tmode: parameter-level Hankinson" << endln;
        s << "\tmaterial: " << theMat->getTag() << endln;
        for (int i = 0; i < numPar; i++)
            s << "\t  " << parNames[i].c_str() << ": " << parVal0[i]
              << " (0 deg) -- " << parVal90[i] << " (90 deg)" << endln;
        s << "\talpha (deg): " << (n1_t != 0.0 || n0_t != 0.0 ?
              atan2(fabs(n1_t), fabs(n0_t)) * 180.0 / 3.14159265358979323846 : 0.0) << endln;
    }
}

// ---------------------------------------------------------------------
// Responses / parameters
// ---------------------------------------------------------------------

Response *
HankinsonSection2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc < 1)
        return SectionForceDeformation::setResponse(argv, argc, output);

    if (strcmp(argv[0], "alpha") == 0)
        return new MaterialResponse(this, 100, 0.0);
    if (strcmp(argv[0], "alphaDeg") == 0)
        return new MaterialResponse(this, 101, 0.0);
    if (strcmp(argv[0], "delta") == 0)
        return new MaterialResponse(this, 102, 0.0);
    if (strcmp(argv[0], "slip") == 0)
        return new MaterialResponse(this, 103, 0.0);
    if (strcmp(argv[0], "sign") == 0)
        return new MaterialResponse(this, 104, 0.0);
    if (strcmp(argv[0], "force") == 0)
        return new MaterialResponse(this, 105, 0.0);
    if (strcmp(argv[0], "stiffness") == 0)
        return new MaterialResponse(this, 106, 0.0);
    if (strcmp(argv[0], "nref") == 0)
        return new MaterialResponse(this, 107, Vector(2));
    if (strcmp(argv[0], "state") == 0)
        return new MaterialResponse(this, 108, Vector(7));
    if (strcmp(argv[0], "explicitError") == 0)
        return new MaterialResponse(this, 109, 0.0);

    if (strcmp(argv[0], "material") == 0 && argc > 1)
        return theMat->setResponse(&argv[1], argc - 1, output);

    return SectionForceDeformation::setResponse(argv, argc, output);
}

int
HankinsonSection2d::getResponse(int responseID, Information &info)
{
    double alpha = atan2(fabs(n1_t), fabs(n0_t));
    if (n0_t == 0.0 && n1_t == 0.0) alpha = 0.0;

    switch (responseID) {
    case 100: return info.setDouble(alpha);
    case 101: return info.setDouble(alpha * 180.0 / 3.14159265358979323846);
    case 102: return info.setDouble(delta_t);
    case 103: return info.setDouble(s_t);
    case 104: return info.setDouble((double)sigma_t);
    case 105: return info.setDouble(F_t);
    case 106: return info.setDouble(k_t);
    case 107: {
        Vector v(2); v(0) = n0_t; v(1) = n1_t;
        return info.setVector(v);
    }
    case 108: {
        Vector v(7);
        v(0) = alpha; v(1) = delta_t; v(2) = s_t; v(3) = (double)sigma_t;
        v(4) = F_t; v(5) = k_t; v(6) = (s_t != 0.0) ? F_t / s_t : 0.0;
        return info.setVector(v);
    }
    case 109: return info.setDouble(CexplicitErr);
    default:
        return SectionForceDeformation::getResponse(responseID, info);
    }
}

int
HankinsonSection2d::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    if (strcmp(argv[0], "material") == 0 && argc > 1)
        return theMat->setParameter(&argv[1], argc - 1, param);

    return -1;
}

int
HankinsonSection2d::updateParameter(int parameterID, Information &info)
{
    return -1;
}
