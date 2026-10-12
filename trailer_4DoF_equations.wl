(* Two-wheeled trailer, 4 DoF model: nonlinear governing equations
     MassMatrix.{psi'', theta'', phi'', u''} + c = 0
   Usage in Mathematica:
     eqs = Get["trailer_4DoF_equations.wl"];
     vals = Join[{v -> 30, Kdpsi -> 5000, psidDelayed -> 0}, eqs["ReferenceParameters"]];
     state = {psi -> 0.01, theta -> 0, phi -> 0, u -> 0, psid -> 0, thetad -> 0, phid -> 0, ud -> 0};
     laws = ReplaceRepeated[ReplaceAll[Join[eqs["ForceLaws"], eqs["Kinematics"]], state], vals];
     Mnum = ReplaceRepeated[ReplaceAll[eqs["MassMatrix"], state], vals];
     cnum = ReplaceRepeated[ReplaceRepeated[ReplaceAll[eqs["c"], state], vals], laws];
     qdd = LinearSolve[Mnum, -cnum]
*)

<|"Description" -> "Two-wheeled trailer, 4 DoF model (H. Z. Horvath, D. Takacs, M. Antali, 2026). \
Equations of motion: MassMatrix.{psi'', theta'', phi'', u''} + c = 0, see Eq. (43) of the paper. \
The forces in c are symbols, defined in ForceLaws in terms of the quantities in \
Kinematics. R: right wheel, L: left wheel. psidDelayed is the delayed yaw rate psid(t - tau) used \
by the braking forces. The parameters v (towing speed) and Kdpsi (control gain) are not included in \
ReferenceParameters.", "GeneralizedCoordinates" -> {psi, theta, phi, u}, 
 "GeneralizedVelocities" -> {psid, thetad, phid, ud}, 
 "MassMatrix" -> 
  {{((JCy + JCz + 2*a^2*m + h^2*m - 4*a*l*m + 2*l^2*m - (JCy - JCz + h^2*m)*Cos[2*phi])*
       Cos[theta]^2)/2 + (JCx + h^2*m)*Sin[theta]^2 + h*(a - l)*m*Cos[phi]*Sin[2*theta], 
    Sin[phi]*((JCy - JCz + h^2*m)*Cos[phi]*Cos[theta] + h*(-a + l)*m*Sin[theta]), 
    -(h*(a - l)*m*Cos[phi]*Cos[theta]) - (JCx + h^2*m)*Sin[theta], 
    m*((a - l)*Cos[psi]*Cos[theta] + h*(Sin[phi]*Sin[psi] + Cos[phi]*Cos[psi]*Sin[theta]))}, 
   {Sin[phi]*((JCy - JCz + h^2*m)*Cos[phi]*Cos[theta] + h*(-a + l)*m*Sin[theta]), 
    (JCy + a^2*m + h^2*m - 2*a*l*m + l^2*m)*Cos[phi]^2 + (JCz + a^2*m - 2*a*l*m + l^2*m)*
      Sin[phi]^2, h*(a - l)*m*Sin[phi], m*Sin[psi]*(h*Cos[phi]*Cos[theta] + (-a + l)*Sin[theta])}, 
   {-(h*(a - l)*m*Cos[phi]*Cos[theta]) - (JCx + h^2*m)*Sin[theta], h*(a - l)*m*Sin[phi], 
    JCx + h^2*m, -(h*m*(Cos[phi]*Cos[psi] + Sin[phi]*Sin[psi]*Sin[theta]))}, 
   {m*((a - l)*Cos[psi]*Cos[theta] + h*(Sin[phi]*Sin[psi] + Cos[phi]*Cos[psi]*Sin[theta])), 
    m*Sin[psi]*(h*Cos[phi]*Cos[theta] + (-a + l)*Sin[theta]), 
    -(h*m*(Cos[phi]*Cos[psi] + Sin[phi]*Sin[psi]*Sin[theta])), m}}, 
 "c" -> {h*(a - l)*m*phid^2*Cos[theta]*Sin[phi] + 2*phid*psid*Cos[theta]*Sin[phi]*
     ((JCy - JCz + h^2*m)*Cos[phi]*Cos[theta] + h*(-a + l)*m*Sin[theta]) - 
    thetad^2*Sin[phi]*(h*(a - l)*m*Cos[theta] + (JCy - JCz + h^2*m)*Cos[phi]*Sin[theta]) + 
    (thetad*(4*phid*(-JCx - h^2*m + (JCy - JCz + h^2*m)*Cos[2*phi])*Cos[theta] + 
       psid*(4*h*(a - l)*m*Cos[phi - 2*theta] + 4*h*(a - l)*m*Cos[phi + 2*theta] - 
         2*(-2*JCx + JCy + JCz + 2*a^2*m - h^2*m - 4*a*l*m + 2*l^2*m - (JCy - JCz + h^2*m)*
            Cos[2*phi])*Sin[2*theta])))/4 + 
    (Sec[phi]*Sec[theta]*(6*FLlat*l + 6*FRlat*l + 2*FLlat*l*Cos[2*phi] + 2*FRlat*l*Cos[2*phi] + 
       2*FLlat*w*Cos[phi - 2*theta] - 2*FRlat*w*Cos[phi - 2*theta] - 
       4*(FLbrake - FRbrake)*w*Cos[phi - theta] - 2*(FLlat + FRlat)*l*Cos[2*theta] + 
       FLlat*l*Cos[2*(-phi + theta)] + FRlat*l*Cos[2*(-phi + theta)] - 
       4*FLbrake*w*Cos[phi + theta] + 4*FRbrake*w*Cos[phi + theta] + FLlat*l*Cos[2*(phi + theta)] + 
       FRlat*l*Cos[2*(phi + theta)] - 2*FLlat*w*Cos[phi + 2*theta] + 2*FRlat*w*Cos[phi + 2*theta] - 
       2*FLbrake*h0*Sin[2*phi] - 2*FRbrake*h0*Sin[2*phi] + 8*FLlat*h0*Sin[theta] + 
       8*FRlat*h0*Sin[theta] + FLbrake*h0*Sin[2*(-phi + theta)] + 
       FRbrake*h0*Sin[2*(-phi + theta)] - FLbrake*h0*Sin[2*(phi + theta)] - 
       FRbrake*h0*Sin[2*(phi + theta)]))/(8*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]), 
   (8*a*h*m*phid^2*Cos[phi] - 8*h*l*m*phid^2*Cos[phi] - 2*a*g*m*Cos[theta] + 2*g*l*m*Cos[theta] + 
     8*JCx*phid*psid*Cos[theta] + 8*h^2*m*phid*psid*Cos[theta] + 8*JCy*phid*psid*Cos[2*phi]*
      Cos[theta] - 8*JCz*phid*psid*Cos[2*phi]*Cos[theta] + 8*h^2*m*phid*psid*Cos[2*phi]*
      Cos[theta] - 8*a*h*m*psid^2*Cos[phi]*Cos[theta]^2 + 8*h*l*m*psid^2*Cos[phi]*Cos[theta]^2 - 
     6*a*g*m*Sec[theta] + 6*g*l*m*Sec[theta] - 8*FLsusp*l*Sec[phi]*Sec[theta]^2 - 
     8*FRsusp*l*Sec[phi]*Sec[theta]^2 - 8*JCy*phid*thetad*Sin[2*phi] + 
     8*JCz*phid*thetad*Sin[2*phi] - 8*h^2*m*phid*thetad*Sin[2*phi] - 3*g*h*m*Cos[phi]*Sin[theta] - 
     16*a*h*m*phid*psid*Cos[phi]*Sin[theta] + 16*h*l*m*phid*psid*Cos[phi]*Sin[theta] - 
     3*g*h*m*Sec[phi]*Sin[theta] + 8*a*h*m*psid^2*Cos[phi]*Sin[theta]^2 - 
     8*h*l*m*psid^2*Cos[phi]*Sin[theta]^2 - (2*FLbrake*h0*Cos[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (2*FRbrake*h0*Cos[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + (2*FLbrake*h0*Sec[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + (2*FRbrake*h0*Sec[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (2*FLbrake*h0*Cos[phi]*Sec[theta]^2)/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (2*FRbrake*h0*Cos[phi]*Sec[theta]^2)/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (6*FLbrake*h0*Sec[phi]*Sec[theta]^2)/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (6*FRbrake*h0*Sec[phi]*Sec[theta]^2)/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (4*FLlat*l*Sin[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (4*FRlat*l*Sin[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (4*FLlat*l*Sec[theta]^2*Sin[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (4*FRlat*l*Sec[theta]^2*Sin[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 4*JCx*psid^2*Sin[2*theta] + 
     2*JCy*psid^2*Sin[2*theta] + 2*JCz*psid^2*Sin[2*theta] + 4*a^2*m*psid^2*Sin[2*theta] - 
     2*h^2*m*psid^2*Sin[2*theta] - 8*a*l*m*psid^2*Sin[2*theta] + 4*l^2*m*psid^2*Sin[2*theta] - 
     2*JCy*psid^2*Cos[phi]^2*Sin[2*theta] + 2*JCz*psid^2*Cos[phi]^2*Sin[2*theta] - 
     2*h^2*m*psid^2*Cos[phi]^2*Sin[2*theta] + 2*JCy*psid^2*Sin[phi]^2*Sin[2*theta] - 
     2*JCz*psid^2*Sin[phi]^2*Sin[2*theta] + 2*h^2*m*psid^2*Sin[phi]^2*Sin[2*theta] + 
     3*g*h*m*Sin[phi]*Sin[theta]*Tan[phi] + (2*FLbrake*h0*Sin[phi]*Tan[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + (2*FRbrake*h0*Sin[phi]*Tan[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + (2*FLbrake*h0*Sec[theta]^2*Sin[phi]*Tan[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + (2*FRbrake*h0*Sec[theta]^2*Sin[phi]*Tan[phi])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - g*h*m*Cos[phi]*Sec[theta]*Tan[theta] - 
     8*FLsusp*h0*Sec[phi]*Sec[theta]*Tan[theta] - 8*FRsusp*h0*Sec[phi]*Sec[theta]*Tan[theta] - 
     g*h*m*Sec[phi]*Sec[theta]*Tan[theta] + 6*a*g*m*Sin[theta]*Tan[theta] - 
     6*g*l*m*Sin[theta]*Tan[theta] - (8*FLbrake*l*Sec[phi]*Sec[theta]*Tan[theta])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - (8*FRbrake*l*Sec[phi]*Sec[theta]*Tan[theta])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + g*h*m*Sec[theta]*Sin[phi]*Tan[phi]*Tan[theta] + 
     g*h*m*Cos[phi]*Sin[theta]*Tan[theta]^2 + g*h*m*Sec[phi]*Sin[theta]*Tan[theta]^2 + 
     (2*FLbrake*h0*Cos[phi]*Tan[theta]^2)/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + 
     (2*FRbrake*h0*Cos[phi]*Tan[theta]^2)/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     (2*FLbrake*h0*Sec[phi]*Tan[theta]^2)/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     (2*FRbrake*h0*Sec[phi]*Tan[theta]^2)/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + 
     (4*FLlat*l*Sin[phi]*Tan[theta]^2)/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + 
     (4*FRlat*l*Sin[phi]*Tan[theta]^2)/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     g*h*m*Sin[phi]*Sin[theta]*Tan[phi]*Tan[theta]^2 - 
     (2*FLbrake*h0*Sin[phi]*Tan[phi]*Tan[theta]^2)/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     (2*FRbrake*h0*Sin[phi]*Tan[phi]*Tan[theta]^2)/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2])/8, 
   (-8*FLlat*h0*Sec[phi]*Sec[theta] - 8*FRlat*h0*Sec[phi]*Sec[theta] - 
     8*FLsusp*w*Sec[phi]^2*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + 
     8*FRsusp*w*Sec[phi]^2*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     3*g*h*m*Cos[theta]*Sin[phi]*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + 
     4*(JCy - JCz + h^2*m)*thetad^2*Sin[2*phi]*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     g*h*m*Sec[phi]^2*Sec[theta]*Sin[3*phi]*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     8*psid^2*Cos[theta]*Sin[phi]*((JCy - JCz + h^2*m)*Cos[phi]*Cos[theta] + 
       h*(-a + l)*m*Sin[theta])*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     4*psid*thetad*(2*(JCx + h^2*m + (JCy - JCz + h^2*m)*Cos[2*phi])*Cos[theta] + 
       4*h*(-a + l)*m*Cos[phi]*Sin[theta])*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] - 
     8*FLlat*w*Tan[phi] + 8*FRlat*w*Tan[phi] + 4*FLbrake*l*Cos[theta]*Sec[phi]*Tan[phi] + 
     4*FRbrake*l*Cos[theta]*Sec[phi]*Tan[phi] - 4*FLbrake*l*Sec[phi]*Sec[theta]*Tan[phi] - 
     4*FRbrake*l*Sec[phi]*Sec[theta]*Tan[phi] - g*h*m*Cos[theta]*Sec[phi]*
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi] - 8*FLsusp*h0*Sec[phi]*Sec[theta]*
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi] - 8*FRsusp*h0*Sec[phi]*Sec[theta]*
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi] - g*h*m*Sec[phi]*Sec[theta]*
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi] - 8*FLbrake*w*Sin[theta]*Tan[phi]^2 + 
     8*FRbrake*w*Sin[theta]*Tan[phi]^2 + g*h*m*Cos[theta]*Sin[phi]*
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi]^2 - 8*FLlat*l*Sec[phi]*Tan[theta] - 
     8*FRlat*l*Sec[phi]*Tan[theta] - 8*FLbrake*h0*Sec[phi]*Tan[phi]*Tan[theta] - 
     8*FRbrake*h0*Sec[phi]*Tan[phi]*Tan[theta] - 4*FLbrake*l*Sec[phi]*Sin[theta]*Tan[phi]*
      Tan[theta] - 4*FRbrake*l*Sec[phi]*Sin[theta]*Tan[phi]*Tan[theta] - 
     8*FLsusp*l*Sec[phi]*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi]*Tan[theta] - 
     8*FRsusp*l*Sec[phi]*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi]*Tan[theta] + 
     2*g*h*m*Sec[phi]*Sin[theta]*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi]*Tan[theta] + 
     2*g*h*m*Cos[2*phi]*Sec[phi]*Sin[theta]*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*Tan[phi]*
      Tan[theta])/(8*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]), 
   -(m*thetad^2*Sin[psi]*((a - l)*Cos[theta] + h*Cos[phi]*Sin[theta])) - 
    2*m*thetad*(h*phid*Cos[theta]*Sin[phi]*Sin[psi] - psid*Cos[psi]*(h*Cos[phi]*Cos[theta] + 
        (-a + l)*Sin[theta])) + (-(FLlat*Cos[phi]*Cos[psi]) - FRlat*Cos[phi]*Cos[psi] + 
      FLbrake*Cos[phi]*Sin[psi] + FRbrake*Cos[phi]*Sin[psi] - FLbrake*Cos[psi]*Sin[phi]*
       Sin[theta] - FRbrake*Cos[psi]*Sin[phi]*Sin[theta] - FLlat*Sin[phi]*Sin[psi]*Sin[theta] - 
      FRlat*Sin[phi]*Sin[psi]*Sin[theta] - FAy0*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + 
      2*h*m*phid*psid*(Cos[phi]*Sin[psi] - Cos[psi]*Sin[phi]*Sin[theta])*
       Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + h*m*phid^2*(Cos[psi]*Sin[phi] - 
        Cos[phi]*Sin[psi]*Sin[theta])*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2] + 
      m*psid^2*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]*(h*Cos[psi]*Sin[phi] - 
        Sin[psi]*((a - l)*Cos[theta] + h*Cos[phi]*Sin[theta])))/
     Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]}, 
 "ForceLaws" -> {FRlat -> muR*NR, FLlat -> muL*NL, 
   NR -> (Sec[phi]*Sec[theta]*(FRbrake*Sin[theta] + 
       FRsusp*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]))/(-(muR*Cos[theta]*Sin[phi]) + 
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]), 
   NL -> (Sec[phi]*Sec[theta]*(FLbrake*Sin[theta] + 
       FLsusp*Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]))/(-(muL*Cos[theta]*Sin[phi]) + 
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2]), 
   muR -> Dmagic*Sin[Cmagic*ArcTan[alphaR*Bmagic - Emagic*(alphaR*Bmagic - 
           ArcTan[alphaR*Bmagic])]], 
   muL -> Dmagic*Sin[Cmagic*ArcTan[alphaL*Bmagic - Emagic*(alphaL*Bmagic - 
           ArcTan[alphaL*Bmagic])]], 
   alphaR -> -ArcTan[vRprime0 . eLat/Sqrt[(vRprime0 - eLat*vRprime0 . eLat) . 
         (vRprime0 - eLat*vRprime0 . eLat)]], 
   alphaL -> -ArcTan[vLprime0 . eLat/Sqrt[(vLprime0 - eLat*vLprime0 . eLat) . 
         (vLprime0 - eLat*vLprime0 . eLat)]], 
   FRsusp -> Max[0, (-(k*(-l0 + sR)) - c*sRdot)*UnitStep[lmax - sR] + 
      (-(kInf*(-lmin + sR)) - cInf*sRdot)*UnitStep[lmin - sR]], 
   FLsusp -> Max[0, (-(k*(-l0 + sL)) - c*sLdot)*UnitStep[lmax - sL] + 
      (-(kInf*(-lmin + sL)) - cInf*sLdot)*UnitStep[lmin - sL]], 
   FRbrake -> Kdpsi*Max[0, -Omegadz + psidDelayed], 
   FLbrake -> Kdpsi*Max[0, -Omegadz - psidDelayed], FAy0 -> -(klat*u) - clat*ud}, 
 "Kinematics" -> {eLat -> {(-(Cos[phi]*Sin[psi]) + Cos[psi]*Sin[phi]*Sin[theta])/
      Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2], 
     (Cos[phi]*Cos[psi] + Sin[phi]*Sin[psi]*Sin[theta])/Sqrt[Cos[phi]^2 + Sin[phi]^2*Sin[theta]^2], 
     0}, vRprime0 -> {v + psid*w*Cos[phi]*Cos[psi] + l*psid*Cos[theta]*Sin[psi] - 
      h0*psid*Cos[psi]*Sec[theta]*Tan[phi] + psid*w*Cos[psi]*Sin[phi]*Tan[phi] - 
      thetad*Sec[theta]^2*(Cos[psi]*(h0 + l*Sin[theta]) + Sin[psi]*(l + h0*Sin[theta])*Tan[phi]) + 
      h0*psid*Sin[psi]*Tan[theta] + l*psid*Sin[psi]*Sin[theta]*Tan[theta] - 
      l*psid*Cos[psi]*Tan[phi]*Tan[theta] - phid*Sec[phi]^2*Sin[psi]*(h0*Sec[theta] - w*Sin[phi] + 
        l*Tan[theta]), ud - l*psid*Cos[psi]*Cos[theta] + h0*phid*Cos[psi]*Sec[theta] - 
      phid*w*Cos[psi]*Sin[phi] + psid*w*Cos[phi]*Sin[psi] - h0*psid*Sec[theta]*Sin[psi]*Tan[phi] + 
      psid*w*Sin[phi]*Sin[psi]*Tan[phi] + h0*phid*Cos[psi]*Sec[theta]*Tan[phi]^2 - 
      phid*w*Cos[psi]*Sin[phi]*Tan[phi]^2 - thetad*Sec[theta]^2*(Sin[psi]*(h0 + l*Sin[theta]) - 
        Cos[psi]*(l + h0*Sin[theta])*Tan[phi]) + l*phid*Cos[psi]*Tan[theta] - 
      h0*psid*Cos[psi]*Tan[theta] - l*psid*Cos[psi]*Sin[theta]*Tan[theta] - 
      l*psid*Sin[psi]*Tan[phi]*Tan[theta] + l*phid*Cos[psi]*Tan[phi]^2*Tan[theta], 0}, 
   vLprime0 -> {v - psid*w*Cos[phi]*Cos[psi] + l*psid*Cos[theta]*Sin[psi] - 
      h0*psid*Cos[psi]*Sec[theta]*Tan[phi] - psid*w*Cos[psi]*Sin[phi]*Tan[phi] - 
      thetad*Sec[theta]^2*(Cos[psi]*(h0 + l*Sin[theta]) + Sin[psi]*(l + h0*Sin[theta])*Tan[phi]) + 
      h0*psid*Sin[psi]*Tan[theta] + l*psid*Sin[psi]*Sin[theta]*Tan[theta] - 
      l*psid*Cos[psi]*Tan[phi]*Tan[theta] - phid*Sec[phi]^2*Sin[psi]*(h0*Sec[theta] + w*Sin[phi] + 
        l*Tan[theta]), ud - l*psid*Cos[psi]*Cos[theta] + h0*phid*Cos[psi]*Sec[theta] + 
      phid*w*Cos[psi]*Sin[phi] - psid*w*Cos[phi]*Sin[psi] - h0*psid*Sec[theta]*Sin[psi]*Tan[phi] - 
      psid*w*Sin[phi]*Sin[psi]*Tan[phi] + h0*phid*Cos[psi]*Sec[theta]*Tan[phi]^2 + 
      phid*w*Cos[psi]*Sin[phi]*Tan[phi]^2 - thetad*Sec[theta]^2*(Sin[psi]*(h0 + l*Sin[theta]) - 
        Cos[psi]*(l + h0*Sin[theta])*Tan[phi]) + l*phid*Cos[psi]*Tan[theta] - 
      h0*psid*Cos[psi]*Tan[theta] - l*psid*Cos[psi]*Sin[theta]*Tan[theta] - 
      l*psid*Sin[psi]*Tan[phi]*Tan[theta] + l*phid*Cos[psi]*Tan[phi]^2*Tan[theta], 0}, 
   sR -> Sec[phi]*(h0*Sec[theta] - w*Sin[phi] + l*Tan[theta]), 
   sL -> Sec[phi]*(h0*Sec[theta] + w*Sin[phi] + l*Tan[theta]), 
   sRdot -> -(phid*w) + l*thetad*Cos[phi] - (thetad*Sec[phi]*Sec[theta]^2*
       (-6*l + 2*l*Cos[2*phi] + l*Cos[2*(phi - theta)] + 2*l*Cos[2*theta] + 
        l*Cos[2*(phi + theta)] - 8*h0*Sin[theta]))/8 + phid*Sec[phi]*Tan[phi]*
      (h0*Sec[theta] - w*Sin[phi] + l*Tan[theta]), 
   sLdot -> phid*w + l*thetad*Cos[phi] - (thetad*Sec[phi]*Sec[theta]^2*
       (-6*l + 2*l*Cos[2*phi] + l*Cos[2*(phi - theta)] + 2*l*Cos[2*theta] + 
        l*Cos[2*(phi + theta)] - 8*h0*Sin[theta]))/8 + phid*Sec[phi]*Tan[phi]*
      (h0*Sec[theta] + w*Sin[phi] + l*Tan[theta])}, "ReferenceParameters" -> 
  {a -> 0.24, h -> 0.21, h0 -> 0.35, l -> 3.77, w -> 0.95, m -> 879, g -> 9.81, JCx -> 554, 
   JCy -> 2107, JCz -> 2601, k -> 30000, c -> 950, kInf -> 2*k, cInf -> 2*c, klat -> 13000, 
   clat -> 750, l0 -> (2*h0*k*l - a*g*m + g*l*m)/(2*k*l), lmin -> 0.28, lmax -> 0.42, 
   Bmagic -> 32.45, Cmagic -> 1.25, Dmagic -> 0.65, Emagic -> -1.24, Omegadz -> 0.1}|>
