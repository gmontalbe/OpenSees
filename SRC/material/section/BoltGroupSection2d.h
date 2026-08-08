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

// Description: This file contains the class definition for BoltGroupSection2d.
//
// BoltGroupSection2d turns a group of mechanical fasteners (dowels, bolts,
// screws, ...) into a SINGLE 3-DOF SectionForceDeformation: it maps the
// (N, V, M) deformation of the connection's reference point RP onto each
// fastener's own 2-D in-plane slip by that fastener's lever arm, evaluates
// one order-2 sub-section per fastener, and sums the results back into the
// (N, V, M) stress resultant and a 3x3 tangent.
//
//     e   = (u_x, u_y, theta)                 at the reference point
//     B_i = [[1, 0, -c*y_i],                  fastener i at offset (x_i, y_i)
//            [0, 1, +c*x_i]]                  in SECTION-LOCAL axes
//     e_i = B_i * e     ->  sub-section i  ->  q_i (2), k_i (2x2)
//     q   = sum_i B_i^T q_i                   (N, V, M) at RP
//     k   = sum_i B_i^T k_i B_i               3x3
//
// This is the same assembly a FiberSection performs, with a 2-DOF coupled
// "fiber" (a whole fastener) in place of a uniaxial one, and it is the
// literal formulation of the reference (monotonic-only) Python
// implementation this generalizes -- see
// connection-mn-interaction/schweigler_implementation, function
// ModelFunctions.internalForces, which sums Fx, Fy and
// M = sum(F_iy*x_i) - sum(F_ix*y_i) over exactly this geometry.
//
// The sub-section is supplied as a PROTOTYPE and copied once per fastener,
// so every fastener carries its own fully independent state and history.
// Any order-2 section works; the intended one is HankinsonSection2d, which
// makes each fastener a Schweigler-style dowel whose scalar slip law has
// its named parameters Hankinson-interpolated between a 0-degree and a
// 90-degree calibration by that fastener's OWN evolving load-to-grain
// angle. Since every fastener in a loaded group sits at a different angle
// alpha, the group's response is not reducible to any single section --
// which is precisely why this class exists.
//
// Relation to the equivalent hand-built model: the same connection can be
// assembled out of one zeroLengthSection element per fastener, all tied to
// a shared master node with rigidLink('beam', ...) -- see
// EXAMPLES/ExamplePython/HankinsonDowelGroup.py. That model and this class
// are two spellings of the same mechanics and (with ROT_LINEAR) agree to
// round-off: rigidLink builds an MP_Constraint with a CONSTANT
// transformation matrix, i.e. it is itself linear in the rotation, exactly
// like B_i above with c == 1. This class replaces N elements plus N
// multi-point constraints with one element and 3 degrees of freedom.
//
// ROTATION KINEMATICS (-rot):
//   ROT_LINEAR (default)  theta is used directly, c == 1, B_i constant.
//   ROT_TAN               tan(theta) is used, as ModelFunctions.dowelDisp
//                         does, and c = sec^2(theta) = 1 + tan^2(theta).
//                         The moment then acquires that same factor, since
//                         the work conjugate of theta under this map is
//                         c * sum_i (x_i q_iy - y_i q_ix), and the tangent
//                         acquires the geometric term 2*tan(theta)*c*m in
//                         k(2,2) -- see setTrialSectionDeformation. Note
//                         that the reference implementation applies the
//                         tan() map to the slips but NOT the matching
//                         factor to its moment, so it is not energy
//                         consistent; ROT_TAN therefore reproduces the
//                         reference's slips exactly and its moment only to
//                         O(theta^2). The two rotation modes likewise agree
//                         to O(theta^2), which at the rotations dowel-type
//                         connections reach is below the staggering error
//                         the sub-section's own lagged Hankinson angle
//                         already introduces (HankinsonSection2d::commitState).
//
// CONTACT ROWS (-contact) are optional compression/bearing springs at
// heights y_j that act along the SECTION-LOCAL X ONLY -- the timber-to-
// timber (or steel-to-timber) contact that closes on one side of a
// connection in bending. Each row is a copy of one UniaxialMaterial (a
// compression-only law such as ElasticPPGap) with its own independent
// history, driven by u_x - c'*y_j, and contributes
//
//     B_j = [1, 0, -c*y_j]     q += B_j^T F_j     k += B_j^T E_j B_j
//
// mirroring ModelFunctions.contactDisp/contactForce. The material is a
// force-vs-displacement law: any tributary area or contact length must
// already be baked into it.
//
// The section reports codes {P, Vy, Mz}, is used with the existing
// ZeroLengthSection element, and is oriented entirely by that element's
// '-orient' flag -- the fastener offsets (x_i, y_i) are always expressed in
// the section-local frame, i.e. the element's local x-y' plane. Note that
// ZeroLengthSection only fills the Mz row of its transformation matrix for
// numDOF == 6 and 12, so a 2D model must use ndf=3 and a 3D model ndf=6
// with both local axes given to '-orient'.

#ifndef BoltGroupSection2d_h
#define BoltGroupSection2d_h

#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <vector>

class UniaxialMaterial;

class BoltGroupSection2d : public SectionForceDeformation
{
  public:
    enum RotRule { ROT_LINEAR = 0, ROT_TAN = 1 };

    // theProto must be an order-2 section; it is copied once per fastener
    // and never itself evaluated. x/y are the fastener offsets from the
    // reference point, in section-local axes. contactProto, if given, is
    // copied once per contact row and driven by the local-x slip at y_j.
    BoltGroupSection2d(int tag, SectionForceDeformation &theProto,
                       int numBolts, const double *x, const double *y,
                       int rotRule = ROT_LINEAR,
                       UniaxialMaterial *contactProto = 0,
                       int numContact = 0, const double *yContact = 0);

    BoltGroupSection2d();          // blank object for FEM_ObjectBroker
    ~BoltGroupSection2d();

    const char *getClassType(void) const { return "BoltGroupSection2d"; }

    int setTrialSectionDeformation(const Vector &e);
    const Vector &getSectionDeformation(void) { return e_t; }
    const Vector &getStressResultant(void)    { return q_t; }
    const Matrix &getSectionTangent(void)     { return kb_t; }
    const Matrix &getInitialTangent(void)     { return K0Mat; }

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    SectionForceDeformation *getCopy(void);
    const ID &getType(void);
    int getOrder(void) const { return 3; }

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &info);
    int setParameter(const char **argv, int argc, Parameter &param);

  protected:

  private:
    // ---- configuration (fixed for the life of the object) -------------
    int rotRule;                    // RotRule
    int numBolts;
    std::vector<double> xOff, yOff;

    SectionForceDeformation *theProto;                  // never evaluated
    std::vector<SectionForceDeformation *> theBolts;    // one copy per fastener

    int numContact;
    std::vector<double> yContact;
    UniaxialMaterial *theContactProto;                  // never evaluated
    std::vector<UniaxialMaterial *> theContact;         // one copy per row

    // ---- trial state --------------------------------------------------
    Vector e_t, q_t;
    Matrix kb_t;
    Vector eBolt;      // scratch order-2 deformation handed to a sub-section
    double tRot_t;     // theta or tan(theta), per rotRule
    double cRot_t;     // 1 or sec^2(theta), per rotRule
    double mSum_t;     // sum_i (x_i q_iy - y_i q_ix), before the c factor

    // ---- committed state ----------------------------------------------
    Vector Ce, Cq;
    Matrix Ckb;

    Matrix K0Mat;      // initial tangent, assembled once
    ID     theCode;    // {P, Vy, Mz}

    // Adds fastener i's contribution B_i^T qi to *q and B_i^T ki B_i to *k;
    // any of qi/ki/q/k may be null. Uses the CURRENT cRot_t. This and
    // assembleContact are the only places the B convention documented at
    // the top of this file is spelled out.
    void assembleBolt(int i, const Vector *qi, const Matrix *ki,
                      Vector *q, Matrix *k);
    void assembleContact(int j, const double *Fj, const double *Ej,
                         Vector *q, Matrix *k);

    void buildInitialTangent(void);
    void allocateBolts(void);     // fills theBolts/theContact with copies
    void freeBolts(void);
};

#endif
