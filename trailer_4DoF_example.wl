(* ::Package:: *)

(* Example: evaluation of the nonlinear governing equations of the 4 DoF trailer model *)
(* M(q).q'' + c = 0, where c also depends on the velocities and on the delayed yaw rate *)

(* 1. Load the equations (save this notebook in the folder of the trailer_4DoF_equations.wl file) *)
ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
eqs = Get["trailer_4DoF_equations.wl"];
Keys[eqs]

(* 2. Parameters: towing speed v, control gain Kdpsi, delayed yaw rate psidDelayed = psid(t - tau),
      and the reference parameters of the paper (rules listed first override the later ones) *)
pars = Join[{v -> 30, Kdpsi -> 5000, psidDelayed -> 0}, eqs["ReferenceParameters"]];

(* 3. State: generalized coordinates and velocities *)
state = {psi -> 0.01, theta -> 0, phi -> 0, u -> 0,
   psid -> 0, thetad -> 0, phid -> 0, ud -> 0};

(* 4. Generalized accelerations {psi'', theta'', phi'', u''} at a given state *)
accelerations[state_, pars_] := Module[{laws, M, cvec},
   (* force laws and kinematic quantities evaluated at the given state *)
   laws = ReplaceRepeated[ReplaceAll[Join[eqs["ForceLaws"], eqs["Kinematics"]], state], pars];
   M = ReplaceRepeated[ReplaceAll[eqs["MassMatrix"], state], pars];
   cvec = ReplaceRepeated[ReplaceRepeated[ReplaceAll[eqs["c"], state], pars], laws];
   LinearSolve[N[M], -N[cvec]]];

accelerations[state, pars]

(* 5. Check: the accelerations are zero in the rectilinear motion *)
Chop[accelerations[Thread[{psi, theta, phi, u, psid, thetad, phid, ud} -> 0], pars]]



