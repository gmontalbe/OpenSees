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

// See BoltGroupSection2d.h for a full description.
//
// Tcl / OpenSeesPy command:
//
//   section BoltGroupSection2d $secTag -section $protoTag \
//       (-bolts $n $x1 $y1 $x2 $y2 ... | -grid $nRows $nCols $a1 $a2) \
//       <-contact $matTag $nRows $y1 $y2 ...> <-rot linear|tan>
//
// e.g. the 3x3 dowel connection of Schweigler et al.'s worked example:
//
//   uniaxialMaterial SAWS 20 19.1 1.91 0.05 24460. 0.0189 -0.05 1.0 0.0 0.8 1.05
//   section HankinsonSection2d 10 -paramHankinson 20 \
//       -param F0 19.1 14.18 -param S0 24460. 19640. -param K1 462. 739.
//   section BoltGroupSection2d 1 -section 10 -grid 3 3 0.070 0.030
//   element zeroLengthSection 1 $iNode $jNode 1 -orient 1 0 0

#include <BoltGroupSection2d.h>

#include <UniaxialMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <Response.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// Beyond this the tan() map is meaningless (and near pi/2, singular). A
// dowel connection reaching even a tenth of this is already far past any
// calibration the sub-section could carry.
#define BGS_MAX_THETA 1.5

// ---------------------------------------------------------------------
// Tcl/Python parser
// ---------------------------------------------------------------------

void *
OPS_BoltGroupSection2d()
{
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: section BoltGroupSection2d tag? -section protoTag? "
               << "<-bolts n? x1? y1? ...  |  -grid nRows? nCols? a1? a2?> "
               << "<-contact matTag? n? y1? ...> <-rot linear|tan>\n";
        return 0;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
        opserr << "WARNING invalid section BoltGroupSection2d tag" << endln;
        return 0;
    }

    SectionForceDeformation *proto = 0;
    UniaxialMaterial *contactProto = 0;
    std::vector<double> xs, ys, yc;
    int rotRule = BoltGroupSection2d::ROT_LINEAR;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *flag = OPS_GetString();

        if (strcmp(flag, "-section") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING -section needs a section tag" << endln;
                return 0;
            }
            int protoTag;
            numData = 1;
            if (OPS_GetIntInput(&numData, &protoTag) != 0) {
                opserr << "WARNING invalid -section tag for section "
                       << "BoltGroupSection2d " << tag << endln;
                return 0;
            }
            proto = OPS_getSectionForceDeformation(protoTag);
            if (proto == 0) {
                opserr << "WARNING section " << protoTag
                       << " not found for section BoltGroupSection2d "
                       << tag << endln;
                return 0;
            }
            // Checked here rather than only in the constructor, which can do
            // nothing but exit(): a mistyped prototype tag is ordinary user
            // error and must not take the interpreter down with it.
            if (proto->getOrder() != 2) {
                opserr << "WARNING section " << protoTag << " has order "
                       << proto->getOrder() << ", section BoltGroupSection2d "
                       << tag << " needs an order-2 prototype (one fastener, "
                       << "codes {P, Vy})" << endln;
                return 0;
            }

        } else if (strcmp(flag, "-bolts") == 0 || strcmp(flag, "-dowels") == 0) {
            int n;
            numData = 1;
            if (OPS_GetIntInput(&numData, &n) != 0 || n < 1) {
                opserr << "WARNING " << flag << " needs a positive fastener "
                       << "count, followed by that many x y pairs" << endln;
                return 0;
            }
            if (OPS_GetNumRemainingInputArgs() < 2 * n) {
                opserr << "WARNING " << flag << " " << n << " needs " << 2 * n
                       << " coordinates, only "
                       << OPS_GetNumRemainingInputArgs() << " arguments left"
                       << endln;
                return 0;
            }
            std::vector<double> v(2 * n);
            numData = 2 * n;
            if (OPS_GetDoubleInput(&numData, &v[0]) != 0) {
                opserr << "WARNING invalid " << flag << " coordinates for "
                       << "section BoltGroupSection2d " << tag << endln;
                return 0;
            }
            xs.clear();
            ys.clear();
            for (int i = 0; i < n; i++) {
                xs.push_back(v[2 * i]);
                ys.push_back(v[2 * i + 1]);
            }

        } else if (strcmp(flag, "-grid") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 4) {
                opserr << "WARNING -grid needs nRows, nCols, a1 and a2" << endln;
                return 0;
            }
            int nrc[2];
            numData = 2;
            if (OPS_GetIntInput(&numData, nrc) != 0 || nrc[0] < 1 || nrc[1] < 1) {
                opserr << "WARNING -grid needs positive nRows and nCols" << endln;
                return 0;
            }
            double a[2];
            numData = 2;
            if (OPS_GetDoubleInput(&numData, a) != 0) {
                opserr << "WARNING invalid -grid spacings a1 a2" << endln;
                return 0;
            }
            // Rectangular grid referred to its own centroid, i.e. the
            // connection's reference point. Row-major (row outer, column
            // inner) to match ModelFunctions.dowelCoordinatesRP in the
            // reference implementation, so fastener indices line up with it.
            // That function's edge distances a3/a4 are deliberately absent:
            // they shift every fastener AND the centroid by the same amount,
            // so they cancel exactly out of centroid-referred coordinates.
            xs.clear();
            ys.clear();
            for (int row = 0; row < nrc[0]; row++) {
                for (int col = 0; col < nrc[1]; col++) {
                    xs.push_back(a[0] * (col - 0.5 * (nrc[1] - 1)));
                    ys.push_back(a[1] * (row - 0.5 * (nrc[0] - 1)));
                }
            }

        } else if (strcmp(flag, "-contact") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 2) {
                opserr << "WARNING -contact needs a material tag, a row count "
                       << "and that many y coordinates" << endln;
                return 0;
            }
            int matTag, n;
            numData = 1;
            if (OPS_GetIntInput(&numData, &matTag) != 0) {
                opserr << "WARNING invalid -contact material tag" << endln;
                return 0;
            }
            contactProto = OPS_getUniaxialMaterial(matTag);
            if (contactProto == 0) {
                opserr << "WARNING uniaxialMaterial " << matTag
                       << " not found for -contact of section "
                       << "BoltGroupSection2d " << tag << endln;
                return 0;
            }
            numData = 1;
            if (OPS_GetIntInput(&numData, &n) != 0 || n < 1) {
                opserr << "WARNING -contact needs a positive row count" << endln;
                return 0;
            }
            if (OPS_GetNumRemainingInputArgs() < n) {
                opserr << "WARNING -contact " << n << " needs " << n
                       << " y coordinates, only "
                       << OPS_GetNumRemainingInputArgs() << " arguments left"
                       << endln;
                return 0;
            }
            yc.resize(n);
            numData = n;
            if (OPS_GetDoubleInput(&numData, &yc[0]) != 0) {
                opserr << "WARNING invalid -contact y coordinates" << endln;
                return 0;
            }

        } else if (strcmp(flag, "-rot") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING -rot needs linear or tan" << endln;
                return 0;
            }
            const char *s = OPS_GetString();
            if (strcmp(s, "linear") == 0)
                rotRule = BoltGroupSection2d::ROT_LINEAR;
            else if (strcmp(s, "tan") == 0)
                rotRule = BoltGroupSection2d::ROT_TAN;
            else {
                opserr << "WARNING invalid -rot \"" << s
                       << "\" (want linear|tan)" << endln;
                return 0;
            }

        } else {
            opserr << "WARNING unknown option \"" << flag
                   << "\" for section BoltGroupSection2d " << tag << endln;
            return 0;
        }
    }

    if (proto == 0) {
        opserr << "WARNING section BoltGroupSection2d " << tag
               << " requires a -section prototype" << endln;
        return 0;
    }
    if (xs.empty()) {
        opserr << "WARNING section BoltGroupSection2d " << tag
               << " requires fastener positions (-bolts or -grid)" << endln;
        return 0;
    }

    return new BoltGroupSection2d(tag, *proto, (int)xs.size(), &xs[0], &ys[0],
                                  rotRule, contactProto, (int)yc.size(),
                                  yc.empty() ? 0 : &yc[0]);
}

// ---------------------------------------------------------------------
// Constructors / destructor
// ---------------------------------------------------------------------

BoltGroupSection2d::BoltGroupSection2d(int tag, SectionForceDeformation &proto,
                                       int nBolts, const double *x, const double *y,
                                       int rot, UniaxialMaterial *contactProto,
                                       int nContact, const double *yc)
: SectionForceDeformation(tag, SEC_TAG_BoltGroupSection2d),
  rotRule(rot), numBolts(nBolts), xOff(nBolts), yOff(nBolts),
  theProto(0),
  numContact(contactProto != 0 ? nContact : 0),
  yContact(contactProto != 0 ? nContact : 0),
  theContactProto(0),
  e_t(3), q_t(3), kb_t(3, 3), eBolt(2),
  tRot_t(0.0), cRot_t(1.0), mSum_t(0.0),
  Ce(3), Cq(3), Ckb(3, 3),
  K0Mat(3, 3), theCode(3)
{
    theCode(0) = SECTION_RESPONSE_P;
    theCode(1) = SECTION_RESPONSE_VY;
    theCode(2) = SECTION_RESPONSE_MZ;

    if (proto.getOrder() != 2) {
        opserr << "BoltGroupSection2d::BoltGroupSection2d -- prototype section "
               << proto.getTag() << " has order " << proto.getOrder()
               << ", must be 2\n";
        exit(-1);
    }

    // The group interprets a fastener's deformation/force component 0 as the
    // section-local x (for timber, parallel-to-grain) slip and component 1 as
    // local y, by position and not by code -- so a prototype carrying other
    // codes still works, it just means something different than intended.
    const ID &protoCode = proto.getType();
    if (protoCode(0) != SECTION_RESPONSE_P || protoCode(1) != SECTION_RESPONSE_VY)
        opserr << "WARNING BoltGroupSection2d::BoltGroupSection2d -- prototype "
               << "section " << proto.getTag() << " does not report codes "
               << "{P, Vy}; its components are used positionally regardless\n";

    theProto = proto.getCopy();
    if (theProto == 0) {
        opserr << "BoltGroupSection2d::BoltGroupSection2d -- failed to get a "
               << "copy of the prototype section\n";
        exit(-1);
    }

    for (int i = 0; i < numBolts; i++) {
        xOff[i] = x[i];
        yOff[i] = y[i];
    }

    if (contactProto != 0 && numContact > 0) {
        theContactProto = contactProto->getCopy();
        if (theContactProto == 0) {
            opserr << "BoltGroupSection2d::BoltGroupSection2d -- failed to get "
                   << "a copy of the contact material\n";
            exit(-1);
        }
        for (int j = 0; j < numContact; j++)
            yContact[j] = yc[j];
    }

    this->allocateBolts();
    this->buildInitialTangent();
    this->revertToStart();
}

BoltGroupSection2d::BoltGroupSection2d()
: SectionForceDeformation(0, SEC_TAG_BoltGroupSection2d),
  rotRule(ROT_LINEAR), numBolts(0),
  theProto(0),
  numContact(0), theContactProto(0),
  e_t(3), q_t(3), kb_t(3, 3), eBolt(2),
  tRot_t(0.0), cRot_t(1.0), mSum_t(0.0),
  Ce(3), Cq(3), Ckb(3, 3),
  K0Mat(3, 3), theCode(3)
{
    theCode(0) = SECTION_RESPONSE_P;
    theCode(1) = SECTION_RESPONSE_VY;
    theCode(2) = SECTION_RESPONSE_MZ;
}

BoltGroupSection2d::~BoltGroupSection2d()
{
    this->freeBolts();
    if (theProto != 0)
        delete theProto;
    if (theContactProto != 0)
        delete theContactProto;
}

void
BoltGroupSection2d::allocateBolts(void)
{
    this->freeBolts();

    theBolts.assign(numBolts, (SectionForceDeformation *)0);
    for (int i = 0; i < numBolts; i++) {
        theBolts[i] = theProto->getCopy();
        if (theBolts[i] == 0) {
            opserr << "BoltGroupSection2d::allocateBolts -- failed to get a "
                   << "copy of the prototype section for fastener " << i << endln;
            exit(-1);
        }
        theBolts[i]->revertToStart();
    }

    theContact.assign(numContact, (UniaxialMaterial *)0);
    for (int j = 0; j < numContact; j++) {
        theContact[j] = theContactProto->getCopy();
        if (theContact[j] == 0) {
            opserr << "BoltGroupSection2d::allocateBolts -- failed to get a "
                   << "copy of the contact material for row " << j << endln;
            exit(-1);
        }
        theContact[j]->revertToStart();
    }
}

void
BoltGroupSection2d::freeBolts(void)
{
    for (size_t i = 0; i < theBolts.size(); i++) {
        if (theBolts[i] != 0)
            delete theBolts[i];
    }
    theBolts.clear();

    for (size_t j = 0; j < theContact.size(); j++) {
        if (theContact[j] != 0)
            delete theContact[j];
    }
    theContact.clear();
}

// ---------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------

void
BoltGroupSection2d::assembleBolt(int i, const Vector *qi, const Matrix *ki,
                                 Vector *q, Matrix *k)
{
    // Columns of B_i (2x3): b_0 = (1, 0), b_1 = (0, 1),
    // b_2 = c * (-y_i, x_i) with c = 1 (ROT_LINEAR) or sec^2(theta) (ROT_TAN).
    const double B[3][2] = { { 1.0,                0.0               },
                             { 0.0,                1.0               },
                             { -cRot_t * yOff[i],  cRot_t * xOff[i]  } };

    if (q != 0 && qi != 0) {
        for (int a = 0; a < 3; a++)
            (*q)(a) += B[a][0] * (*qi)(0) + B[a][1] * (*qi)(1);
    }

    if (k != 0 && ki != 0) {
        // k_ab += b_a^T * k_i * b_b, written out in full rather than
        // exploiting symmetry so that a prototype with a non-symmetric
        // tangent assembles correctly too.
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                double s = 0.0;
                for (int p = 0; p < 2; p++)
                    for (int r = 0; r < 2; r++)
                        s += B[a][p] * (*ki)(p, r) * B[b][r];
                (*k)(a, b) += s;
            }
        }
    }
}

void
BoltGroupSection2d::assembleContact(int j, const double *Fj, const double *Ej,
                                    Vector *q, Matrix *k)
{
    // A contact row acts along local x only: B_j = [1, 0, -c*y_j].
    const double B[3] = { 1.0, 0.0, -cRot_t * yContact[j] };

    if (q != 0 && Fj != 0) {
        for (int a = 0; a < 3; a++)
            (*q)(a) += B[a] * (*Fj);
    }

    if (k != 0 && Ej != 0) {
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                (*k)(a, b) += B[a] * (*Ej) * B[b];
    }
}

void
BoltGroupSection2d::buildInitialTangent(void)
{
    // theta == 0 in the undeformed state, so c == 1 in both rotation modes
    // and the geometric term (which carries a factor tan(theta)) vanishes.
    double cSave = cRot_t;
    cRot_t = 1.0;

    K0Mat.Zero();
    for (int i = 0; i < numBolts; i++)
        this->assembleBolt(i, 0, &theBolts[i]->getInitialTangent(), 0, &K0Mat);

    for (int j = 0; j < numContact; j++) {
        double Ej = theContact[j]->getInitialTangent();
        this->assembleContact(j, 0, &Ej, 0, &K0Mat);
    }

    cRot_t = cSave;

    // sum_i k_i * r_i^2 == 0 means every fastener sits at the reference
    // point: the group carries N and V but has no rotational stiffness at
    // all, and the Mz row/column of the tangent is identically zero. That is
    // a legal configuration, but it makes the assembled system singular
    // unless the host element's relative rotation is restrained elsewhere.
    if (K0Mat(2, 2) <= 0.0)
        opserr << "WARNING BoltGroupSection2d -- section " << this->getTag()
               << " has zero rotational stiffness (all " << numBolts
               << " fasteners lie at the reference point); restrain the "
               << "relative rotation of the host element's nodes\n";
}

int
BoltGroupSection2d::setTrialSectionDeformation(const Vector &e)
{
    e_t(0) = e(0);
    e_t(1) = e(1);
    e_t(2) = e(2);

    if (rotRule == ROT_TAN) {
        if (fabs(e(2)) > BGS_MAX_THETA) {
            opserr << "WARNING BoltGroupSection2d -- section " << this->getTag()
                   << ": |theta| = " << fabs(e(2)) << " exceeds " << BGS_MAX_THETA
                   << " rad, where the -rot tan map breaks down\n";
            return -1;
        }
        tRot_t = tan(e(2));
        cRot_t = 1.0 + tRot_t * tRot_t;   // sec^2(theta)
    } else {
        tRot_t = e(2);
        cRot_t = 1.0;
    }

    q_t.Zero();
    kb_t.Zero();
    mSum_t = 0.0;

    int err = 0;
    for (int i = 0; i < numBolts; i++) {
        // e_i = B_i * e -- a rigid-body motion of the reference point,
        // distributed to fastener i by its lever arm.
        eBolt(0) = e(0) - tRot_t * yOff[i];
        eBolt(1) = e(1) + tRot_t * xOff[i];

        err += theBolts[i]->setTrialSectionDeformation(eBolt);

        const Vector &qi = theBolts[i]->getStressResultant();
        this->assembleBolt(i, &qi, &theBolts[i]->getSectionTangent(), &q_t, &kb_t);
        mSum_t += xOff[i] * qi(1) - yOff[i] * qi(0);
    }

    for (int j = 0; j < numContact; j++) {
        err += theContact[j]->setTrialStrain(e(0) - tRot_t * yContact[j]);

        double Fj = theContact[j]->getStress();
        double Ej = theContact[j]->getTangent();
        this->assembleContact(j, &Fj, &Ej, &q_t, &kb_t);
        mSum_t += -yContact[j] * Fj;
    }

    // Geometric term: with tan(theta) kinematics the moment conjugate is
    // c*m, so d(q_M)/d(theta) picks up (dc/dtheta)*m = 2*tan(theta)*c*m on
    // top of the material term already assembled above. It is diagonal, so
    // the tangent stays symmetric.
    if (rotRule == ROT_TAN)
        kb_t(2, 2) += 2.0 * tRot_t * cRot_t * mSum_t;

    return err;
}

// ---------------------------------------------------------------------
// State management
// ---------------------------------------------------------------------

int
BoltGroupSection2d::commitState(void)
{
    int err = 0;
    for (int i = 0; i < numBolts; i++)
        err += theBolts[i]->commitState();
    for (int j = 0; j < numContact; j++)
        err += theContact[j]->commitState();

    if (err != 0)
        return err;

    Ce = e_t;
    Cq = q_t;
    Ckb = kb_t;

    return 0;
}

int
BoltGroupSection2d::revertToLastCommit(void)
{
    int err = 0;
    for (int i = 0; i < numBolts; i++)
        err += theBolts[i]->revertToLastCommit();
    for (int j = 0; j < numContact; j++)
        err += theContact[j]->revertToLastCommit();

    e_t = Ce;
    q_t = Cq;
    kb_t = Ckb;

    if (rotRule == ROT_TAN) {
        tRot_t = tan(Ce(2));
        cRot_t = 1.0 + tRot_t * tRot_t;
    } else {
        tRot_t = Ce(2);
        cRot_t = 1.0;
    }

    return err;
}

int
BoltGroupSection2d::revertToStart(void)
{
    int err = 0;
    for (int i = 0; i < numBolts; i++)
        err += theBolts[i]->revertToStart();
    for (int j = 0; j < numContact; j++)
        err += theContact[j]->revertToStart();

    Ce.Zero();
    Cq.Zero();
    Ckb = K0Mat;

    e_t.Zero();
    q_t.Zero();
    kb_t = K0Mat;

    tRot_t = 0.0;
    cRot_t = 1.0;
    mSum_t = 0.0;

    return err;
}

// ---------------------------------------------------------------------
// Copy / type / serialization
// ---------------------------------------------------------------------

SectionForceDeformation *
BoltGroupSection2d::getCopy(void)
{
    BoltGroupSection2d *theCopy =
        new BoltGroupSection2d(this->getTag(), *theProto, numBolts,
                               numBolts > 0 ? &xOff[0] : 0,
                               numBolts > 0 ? &yOff[0] : 0,
                               rotRule, theContactProto, numContact,
                               numContact > 0 ? &yContact[0] : 0);

    // The constructor above gave the copy a set of VIRGIN fasteners; swap in
    // copies of our own, which carry their committed state. Going through
    // each sub-section's own getCopy() is what keeps the copy independent:
    // for HankinsonSection2d that is what re-binds fresh Parameters onto a
    // fresh material copy, without which every element sharing this section
    // tag would mutate one shared law (see HankinsonSection2d::getCopy).
    for (int i = 0; i < numBolts; i++) {
        SectionForceDeformation *boltCopy = theBolts[i]->getCopy();
        if (boltCopy == 0) {
            opserr << "BoltGroupSection2d::getCopy -- failed to copy fastener "
                   << i << endln;
            delete theCopy;
            return 0;
        }
        delete theCopy->theBolts[i];
        theCopy->theBolts[i] = boltCopy;
    }
    for (int j = 0; j < numContact; j++) {
        UniaxialMaterial *matCopy = theContact[j]->getCopy();
        if (matCopy == 0) {
            opserr << "BoltGroupSection2d::getCopy -- failed to copy contact "
                   << "row " << j << endln;
            delete theCopy;
            return 0;
        }
        delete theCopy->theContact[j];
        theCopy->theContact[j] = matCopy;
    }

    theCopy->Ce = Ce;
    theCopy->Cq = Cq;
    theCopy->Ckb = Ckb;
    theCopy->revertToLastCommit();

    return theCopy;
}

const ID &
BoltGroupSection2d::getType(void)
{
    return theCode;
}

int
BoltGroupSection2d::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    int dbTag = this->getDbTag();

    ID idData(4);
    idData(0) = this->getTag();
    idData(1) = numBolts;
    idData(2) = numContact;
    idData(3) = rotRule;
    res = theChannel.sendID(dbTag, cTag, idData);
    if (res < 0) {
        opserr << "BoltGroupSection2d::sendSelf -- failed to send ID data\n";
        return res;
    }

    // class/db tag pair for the prototype, then for every fastener, then
    // (if any) the contact prototype and every contact row
    int numSec = numBolts + 1;
    int numMat = (numContact > 0) ? numContact + 1 : 0;
    ID classTags(2 * (numSec + numMat));
    for (int i = 0; i < numSec; i++) {
        SectionForceDeformation *theSec = (i == 0) ? theProto : theBolts[i - 1];
        classTags(2 * i) = theSec->getClassTag();
        if (theSec->getDbTag() == 0) {
            int secDbTag = theChannel.getDbTag();
            if (secDbTag != 0)
                theSec->setDbTag(secDbTag);
        }
        classTags(2 * i + 1) = theSec->getDbTag();
    }
    for (int j = 0; j < numMat; j++) {
        UniaxialMaterial *theMat = (j == 0) ? theContactProto : theContact[j - 1];
        int loc = 2 * (numSec + j);
        classTags(loc) = theMat->getClassTag();
        if (theMat->getDbTag() == 0) {
            int matDbTag = theChannel.getDbTag();
            if (matDbTag != 0)
                theMat->setDbTag(matDbTag);
        }
        classTags(loc + 1) = theMat->getDbTag();
    }
    res = theChannel.sendID(dbTag, cTag, classTags);
    if (res < 0) {
        opserr << "BoltGroupSection2d::sendSelf -- failed to send object tags\n";
        return res;
    }

    // fastener offsets, contact offsets, then the committed state
    Vector data(2 * numBolts + numContact + 15);
    for (int i = 0; i < numBolts; i++) {
        data(2 * i)     = xOff[i];
        data(2 * i + 1) = yOff[i];
    }
    int loc = 2 * numBolts;
    for (int j = 0; j < numContact; j++)
        data(loc + j) = yContact[j];
    loc += numContact;
    for (int i = 0; i < 3; i++) {
        data(loc + i)     = Ce(i);
        data(loc + 3 + i) = Cq(i);
    }
    loc += 6;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            data(loc + 3 * i + j) = Ckb(i, j);

    res = theChannel.sendVector(dbTag, cTag, data);
    if (res < 0) {
        opserr << "BoltGroupSection2d::sendSelf -- failed to send Vector data\n";
        return res;
    }

    res = theProto->sendSelf(cTag, theChannel);
    if (res < 0) {
        opserr << "BoltGroupSection2d::sendSelf -- failed to send prototype\n";
        return res;
    }
    for (int i = 0; i < numBolts; i++) {
        res = theBolts[i]->sendSelf(cTag, theChannel);
        if (res < 0) {
            opserr << "BoltGroupSection2d::sendSelf -- failed to send fastener "
                   << i << endln;
            return res;
        }
    }
    if (numContact > 0) {
        res = theContactProto->sendSelf(cTag, theChannel);
        if (res < 0) {
            opserr << "BoltGroupSection2d::sendSelf -- failed to send contact "
                   << "prototype\n";
            return res;
        }
        for (int j = 0; j < numContact; j++) {
            res = theContact[j]->sendSelf(cTag, theChannel);
            if (res < 0) {
                opserr << "BoltGroupSection2d::sendSelf -- failed to send "
                       << "contact row " << j << endln;
                return res;
            }
        }
    }

    return 0;
}

int
BoltGroupSection2d::recvSelf(int cTag, Channel &theChannel,
                             FEM_ObjectBroker &theBroker)
{
    int res = 0;
    int dbTag = this->getDbTag();

    ID idData(4);
    res = theChannel.recvID(dbTag, cTag, idData);
    if (res < 0) {
        opserr << "BoltGroupSection2d::recvSelf -- failed to receive ID data\n";
        return res;
    }
    this->setTag(idData(0));
    numBolts = idData(1);
    numContact = idData(2);
    rotRule = idData(3);

    int numSec = numBolts + 1;
    int numMat = (numContact > 0) ? numContact + 1 : 0;
    ID classTags(2 * (numSec + numMat));
    res = theChannel.recvID(dbTag, cTag, classTags);
    if (res < 0) {
        opserr << "BoltGroupSection2d::recvSelf -- failed to receive object tags\n";
        return res;
    }

    Vector data(2 * numBolts + numContact + 15);
    res = theChannel.recvVector(dbTag, cTag, data);
    if (res < 0) {
        opserr << "BoltGroupSection2d::recvSelf -- failed to receive Vector data\n";
        return res;
    }

    xOff.resize(numBolts);
    yOff.resize(numBolts);
    for (int i = 0; i < numBolts; i++) {
        xOff[i] = data(2 * i);
        yOff[i] = data(2 * i + 1);
    }
    int loc = 2 * numBolts;
    yContact.resize(numContact);
    for (int j = 0; j < numContact; j++)
        yContact[j] = data(loc + j);
    loc += numContact;
    for (int i = 0; i < 3; i++) {
        Ce(i) = data(loc + i);
        Cq(i) = data(loc + 3 + i);
    }
    loc += 6;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Ckb(i, j) = data(loc + 3 * i + j);

    this->freeBolts();
    if (theProto != 0) {
        delete theProto;
        theProto = 0;
    }
    if (theContactProto != 0) {
        delete theContactProto;
        theContactProto = 0;
    }

    theBolts.assign(numBolts, (SectionForceDeformation *)0);
    for (int i = 0; i < numSec; i++) {
        SectionForceDeformation *theSec = theBroker.getNewSection(classTags(2 * i));
        if (theSec == 0) {
            opserr << "BoltGroupSection2d::recvSelf -- could not get a section "
                   << "of class " << classTags(2 * i) << endln;
            return -1;
        }
        theSec->setDbTag(classTags(2 * i + 1));
        res = theSec->recvSelf(cTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "BoltGroupSection2d::recvSelf -- section failed to recvSelf\n";
            delete theSec;
            return res;
        }
        if (i == 0)
            theProto = theSec;
        else
            theBolts[i - 1] = theSec;
    }

    theContact.assign(numContact, (UniaxialMaterial *)0);
    for (int j = 0; j < numMat; j++) {
        int idx = 2 * (numSec + j);
        UniaxialMaterial *theMat = theBroker.getNewUniaxialMaterial(classTags(idx));
        if (theMat == 0) {
            opserr << "BoltGroupSection2d::recvSelf -- could not get a material "
                   << "of class " << classTags(idx) << endln;
            return -1;
        }
        theMat->setDbTag(classTags(idx + 1));
        res = theMat->recvSelf(cTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "BoltGroupSection2d::recvSelf -- material failed to recvSelf\n";
            delete theMat;
            return res;
        }
        if (j == 0)
            theContactProto = theMat;
        else
            theContact[j - 1] = theMat;
    }

    this->buildInitialTangent();
    this->revertToLastCommit();

    return 0;
}

void
BoltGroupSection2d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "BoltGroupSection2d, tag: " << this->getTag() << endln;
        s << "\trotation: " << (rotRule == ROT_TAN ? "tan" : "linear") << endln;
        s << "\tfasteners: " << numBolts << endln;
        for (int i = 0; i < numBolts; i++)
            s << "\t  " << i << ": (" << xOff[i] << ", " << yOff[i] << ")" << endln;
        if (numContact > 0) {
            s << "\tcontact rows: " << numContact << endln;
            for (int j = 0; j < numContact; j++)
                s << "\t  " << j << ": y = " << yContact[j] << endln;
            s << "\tcontact material: " << theContactProto->getClassType()
              << ", tag " << theContactProto->getTag() << endln;
        }
        if (theProto != 0) {
            s << "\tprototype section: " << theProto->getClassType()
              << ", tag " << theProto->getTag() << endln;
            theProto->Print(s, flag);
        }
    }
}

// ---------------------------------------------------------------------
// Responses / parameters
// ---------------------------------------------------------------------

Response *
BoltGroupSection2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc < 1)
        return SectionForceDeformation::setResponse(argv, argc, output);

    // "bolt i <whatever the sub-section understands>" -- the group's
    // replacement for reading one dowel's own zeroLengthSection in the
    // hand-built (one element per fastener) model.
    if (argc > 2 && (strcmp(argv[0], "bolt") == 0 || strcmp(argv[0], "dowel") == 0)) {
        int i = atoi(argv[1]);
        if (i < 0 || i >= numBolts) {
            opserr << "WARNING BoltGroupSection2d::setResponse -- fastener index "
                   << i << " out of range [0, " << numBolts - 1 << "]\n";
            return 0;
        }
        return theBolts[i]->setResponse(&argv[2], argc - 2, output);
    }

    if (argc > 2 && strcmp(argv[0], "contact") == 0) {
        int j = atoi(argv[1]);
        if (j < 0 || j >= numContact) {
            opserr << "WARNING BoltGroupSection2d::setResponse -- contact row "
                   << j << " out of range [0, " << numContact - 1 << "]\n";
            return 0;
        }
        return theContact[j]->setResponse(&argv[2], argc - 2, output);
    }

    // Whole-group vectors, so one recorder covers every fastener.
    if (strcmp(argv[0], "boltForces") == 0 || strcmp(argv[0], "dowelForces") == 0)
        return new MaterialResponse(this, 100, Vector(2 * numBolts));
    if (strcmp(argv[0], "boltSlips") == 0 || strcmp(argv[0], "dowelSlips") == 0)
        return new MaterialResponse(this, 101, Vector(2 * numBolts));
    if (strcmp(argv[0], "centerOfRotation") == 0)
        return new MaterialResponse(this, 102, Vector(2));
    if (numContact > 0 && strcmp(argv[0], "contactForces") == 0)
        return new MaterialResponse(this, 103, Vector(numContact));

    return SectionForceDeformation::setResponse(argv, argc, output);
}

int
BoltGroupSection2d::getResponse(int responseID, Information &info)
{
    switch (responseID) {

    case 100: {
        Vector v(2 * numBolts);
        for (int i = 0; i < numBolts; i++) {
            const Vector &qi = theBolts[i]->getStressResultant();
            v(2 * i)     = qi(0);
            v(2 * i + 1) = qi(1);
        }
        return info.setVector(v);
    }

    case 101: {
        Vector v(2 * numBolts);
        for (int i = 0; i < numBolts; i++) {
            const Vector &ei = theBolts[i]->getSectionDeformation();
            v(2 * i)     = ei(0);
            v(2 * i + 1) = ei(1);
        }
        return info.setVector(v);
    }

    case 102: {
        // The centre of rotation is where the rigid-body slip field vanishes:
        // u_x - t*y = 0 and u_y + t*x = 0, so (x, y) = (-u_y/t, u_x/t). It is
        // undefined for a pure translation (t == 0), reported as (0, 0) --
        // the honest answer there is "at infinity". Compare
        // ModelFunctions.findCR, which locates the same point by intersecting
        // two fasteners' force normals.
        Vector v(2);
        if (tRot_t != 0.0) {
            v(0) = -e_t(1) / tRot_t;
            v(1) =  e_t(0) / tRot_t;
        } else {
            v(0) = 0.0;
            v(1) = 0.0;
        }
        return info.setVector(v);
    }

    case 103: {
        Vector v(numContact);
        for (int j = 0; j < numContact; j++)
            v(j) = theContact[j]->getStress();
        return info.setVector(v);
    }

    default:
        return SectionForceDeformation::getResponse(responseID, info);
    }
}

int
BoltGroupSection2d::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    if (argc > 2 && (strcmp(argv[0], "bolt") == 0 || strcmp(argv[0], "dowel") == 0)) {
        int i = atoi(argv[1]);
        if (i < 0 || i >= numBolts) {
            opserr << "WARNING BoltGroupSection2d::setParameter -- fastener index "
                   << i << " out of range [0, " << numBolts - 1 << "]\n";
            return -1;
        }
        return theBolts[i]->setParameter(&argv[2], argc - 2, param);
    }

    if (argc > 2 && strcmp(argv[0], "contact") == 0) {
        int j = atoi(argv[1]);
        if (j < 0 || j >= numContact) {
            opserr << "WARNING BoltGroupSection2d::setParameter -- contact row "
                   << j << " out of range [0, " << numContact - 1 << "]\n";
            return -1;
        }
        return theContact[j]->setParameter(&argv[2], argc - 2, param);
    }

    // Anything else is offered to every fastener, so that one Parameter can
    // drive the whole group at once (Parameter holds a list of objects).
    int result = -1;
    for (int i = 0; i < numBolts; i++) {
        int ok = theBolts[i]->setParameter(argv, argc, param);
        if (ok >= 0)
            result = ok;
    }

    return result;
}
