(* ::Package:: *)

(* Mathematica Package  *)

(* :Title: CompoundMatrixMethod *)
(* :Author: Simon Pearce <simon.pearce@manchester.ac.uk> *)
(* :Context: CompoundMatrixMethod` *)
(* :Version: 0.7 *)
(* :Date: 2018-05-14 *)

(* :Mathematica Version: 9+ *)
(* :Copyright: (c) 2017-18 Simon Pearce *)

BeginPackage["CompoundMatrixMethod`"];

Unprotect["CompoundMatrixMethod`*"];

Evans::usage = "\
Evans[{k,k0},sys] evaluates the Evans function from the Compound Matrix Method with potential eigenvalue k=k0, for the system defined from ToLinearMatrixForm.
Evans[{k,k0},A,B,C,{x,x0,x1}] evaluates the Evans function from the Compound Matrix Method with k=k0. Here the linear matrix ODE is given by dy/dx=A.y, with B.y=0 at x=x0 and C.y=0 at x=x1.
For either case a complex number is returned, zeroes of this function correspond to zeroes of the original eigenvalue equation.
Evans[{k,k0},{A1,A2},B,C,{F,G},{x,x0,xmatch,x1}] evaluates the Evans function for potential eigenvalue k=k0 for an interface problem.
Here there are two linear matrix ODE is given on the left by dy/dx=A1.y, B.y=0 at x=x0, dz/dx=A2.z, and C.z=0 at x=x1, and the interface conditions are given by
F.y + G.z = 0 at the interface given by x=xmatch.";

ToLinearMatrixForm::usage = "\
ToLinearMatrixForm[eqn,{},depvars,x] takes a list of differential equations in the dependent variables depvars and independent variable x and puts the equations into linear matrix form dy/dx = A.y
ToLinearMatrixForm[eqn,bcs,depvars,{x,x0,x1}] also includes the boundary conditions evaluated at x=x0 and x=x1.
ToLinearMatrixForm[{eqns1,eqns2},{bcs,interfaceEqns},depvars,{x,x0,xmatch,x1}] sets up the function with an interface at x=xmatch, with eqns1 in x0<x<xmatch and eqns2 in xmatch<x<x1.";

Evans::boundarySolutionFailed = "Applying the boundary conditions at `1` failed, perhaps you don't have any conditions there.";
Evans::nonSquareMatrix = "The matrix A is not square.";
Evans::matrixRank = "The matrix A has rank less than its dimensions.";

Evans::boundaryConditionRank = "The rank of the boundary condition matrices does not sum to the size of A";
Evans::nonZeroBoundaryConditions = "Non-zero boundary condition given, ";

Evans::nonNumericalMatrix = "Matrix `1` is not numerical at x = `2` ";
Evans::noEigenvalue = "Warning: Matrix `1` does not contain the potential eigenvalue `2` ";

Evans::incorrectRange = "Range of integration, [`1` `2`] is not numerical with `1`<`2`";

Evans::matchPointIncorrect = "Range of integration, [`1` `2`] does not include the matching point `3`";

Evans::matrixSizesDiffer = "The two matrices `1` and `2` have different sizes.";

Evans::InterfaceTooBig = "Interface problems are only currently supported up to dimension 10, your matrix size is `1`. Please raise an issue on github to ask for larger sizes to be implemented.";

Evans::normalizationConstants = "Normalization constants must be a list of two elements.";

ToLinearMatrixForm::somethingWrong = "Something has gone wrong!";

ToLinearMatrixForm::matrixOnly = "Incorrect number of boundary conditions given (`1` compared to matrix dimension `2`), returning the matrix only";
ToLinearMatrixForm::inhomogeneous = "Inhomogeneous equation, does not work with the Compound Matrix Method";
ToLinearMatrixForm::inhomogeneousBCs = "Boundary values are inhomogeneous, does not work with the Compound Matrix Method (returning vectors as well as matrices) `1` `2`";

ToLinearMatrixForm::linearised = "Original Equation not linear, linearized equation given by `1`";
ToLinearMatrixForm::noLimits = "Please supply limits for the independent variable when providing boundary conditions";

SelectPositiveEigenvectors::usage = "\
SelectPositiveEigenvectors[matrix, x] selects the eigenvectors which correspond to negative real part as x->-inf";

SelectNegativeEigenvectors::usage = "\
SelectNegativeEigenvectors[matrix] selects the eigenvectors which correspond to positive real part at  x->inf";

SelectPositiveEigenvectors::nonNumericalMatrix = "Matrix `1` is not numerical";
SelectNegativeEigenvectors::nonNumericalMatrix = "Matrix `1` is not numerical";

SelectPositiveEigenvectors::noPositiveEigenvectors = "No positive Eigenvectors found";
SelectNegativeEigenvectors::noNegativeEigenvectors = "No negative Eigenvectors found";

CompoundMatrixMethod::usage = "Alternative alias for Evans";
NormalizationConstants::usage = "Normalization constants for the Evans function, for the integration from left and right respectively. Usual settings are either 0 or 1 for both.";

Begin["`Private`"]; (* Begin Private Context *)

(* Initially used CompoundMatrixMethod as the function name, so this is for backward compatibility. *)
CompoundMatrixMethod = Evans;

(* Replacement rule to sort the lists of indices *)
reprules = \[FormalPhi][a_?ListQ] :> Signature[a] \[FormalPhi][Sort[a]];
reprules2= {\[FormalPhi]L[a_?ListQ][q_] :>	Signature[a] \[FormalPhi]L[Sort[a]][q], \[FormalPhi]R[a_?ListQ][q_] :> Signature[a] \[FormalPhi]R[Sort[a]][q]};

(* Generation of the derivatives of the matrix minors, looping over  rule to sort the lists of indices *)
minorsDerivs[list_?VectorQ, len_?NumericQ] := minorsDerivs[list, len] =
    Sum[ Sum[\[FormalCapitalA][y, z] \[FormalPhi][list /. y -> z], {z, Union[Complement[Range[len], list], {y}]}],
	    {y, list}] /. reprules

qMatrix[len_?NumericQ, len2_?NumericQ] := qMatrix[len, len2] =
    Module[{minorsTab},
  minorsTab = Table[minorsDerivs[ii, len], {ii, Subsets[Range[len], {len2}]}]
		  /. Thread[Subsets[Range[len], {len2}] -> Range[Length[Subsets[Range[len], {len2}]]]];
  Transpose[Table[Coefficient[minorsTab, \[FormalPhi][i]], {i, 1, Binomial[len, len2]}]]
]

rulesFG[len_?NumericQ] := rulesFG[len] = {\[FormalPhi]\[FormalPhi]L[{l_}][x_] :>	Sum[F[l, a] \[FormalPhi]L[a][x], {a, 1,	len}],
	\[FormalPhi]\[FormalPhi]L[{l_, m_}][x_] :>	Sum[F[l, a] F[m, b] \[FormalPhi]L[{a, b}][x], {a, 1, len}, {b, 1,	len}],
	\[FormalPhi]\[FormalPhi]L[{l_, m_, n_}][x_] :>	Sum[F[l, a] F[m, b] F[n, c] \[FormalPhi]L[{a, b, c}][x], {a, 1,	len}, {b, 1, len}, {c, 1,	len}],
	\[FormalPhi]\[FormalPhi]L[{l_, m_, n_, o_}][x_] :>	Sum[F[l, a] F[m, b] F[n, c] F[o, d]
		\[FormalPhi]L[{a, b, c, d}][x], {a, 1, len}, {b, 1, len}, {c, 1, len}, {d, 1,	len}],
	\[FormalPhi]\[FormalPhi]L[{l_, m_, n_, o_, p_}][x_] :>	Sum[F[l, a] F[m, b] F[n, c] F[o, d] F[p,e]
		\[FormalPhi]L[{a, b, c, d, e}][x], {a, 1, len}, {b, 1, len}, {c, 1, len}, {d, 1, len}, {e, 1, len}],
	\[FormalPhi]\[FormalPhi]R[{l_}][x_] :>	Sum[G[l, a] \[FormalPhi]R[a][x], {a, 1,	len}],
	\[FormalPhi]\[FormalPhi]R[{l_, m_}][x_] :>	Sum[G[l, a] G[m, b] \[FormalPhi]R[{a, b}][x], {a, 1, len}, {b, 1,	len}],
	\[FormalPhi]\[FormalPhi]R[{l_, m_, n_}][x_] :>	Sum[G[l, a] G[m, b] G[n, c] \[FormalPhi]R[{a, b, c}][x], {a, 1,	len}, {b, 1, len}, {c, 1, len}],
	\[FormalPhi]\[FormalPhi]R[{l_, m_, n_, o_}][x_] :>	Sum[G[l, a] G[m, b] G[n, c] G[o, d] \[FormalPhi]R[{a, b, c, d}][
		x], {a, 1, len}, {b, 1, len}, {c, 1, len}, {d, 1,	len}],
	\[FormalPhi]\[FormalPhi]R[{l_, m_, n_, o_, p_}][x_] :>	Sum[G[l, a] G[m, b] G[n, c] G[o, d] G[p, e]
		\[FormalPhi]R[{a, b, c, d, e}][x], {a, 1, len}, {b, 1, len}, {c, 1, len}, {d, 1, len}, {e, 1, len}]};

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, {AMatrix_?MatrixQ, {}, rightBCMatrix_?MatrixQ, {x_ /; ! NumericQ[x], xa_, xb_,	xm_ : False}},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
		Evans[{\[FormalLambda], \[FormalLambda]0},	AMatrix, N@SelectNegativeEigenvectors[AMatrix /. \[FormalLambda] -> \[FormalLambda]0, x], rightBCMatrix, {x, xa, xb, xm},opts]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, {AMatrix_?MatrixQ, {}, {}, {x_ /; ! NumericQ[x], xa_, xb_, xm_ : False}},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
		Evans[{\[FormalLambda], \[FormalLambda]0}, AMatrix, N@SelectNegativeEigenvectors[AMatrix /. \[FormalLambda] -> \[FormalLambda]0, x],
			N@SelectPositiveEigenvectors[AMatrix /. \[FormalLambda] -> \[FormalLambda]0, x], {x, xa,	xb, xm},opts]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, {AMatrix_?MatrixQ,	leftBCMatrix_?MatrixQ, {}, {x_ /; ! NumericQ[x], xa_, xb_, xm_ : False}},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
		Evans[{\[FormalLambda], \[FormalLambda]0},	AMatrix, leftBCMatrix, N@SelectPositiveEigenvectors[AMatrix /. \[FormalLambda] -> \[FormalLambda]0, x], {x, xa,	xb, xm},opts]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ},	AMatrix_?MatrixQ, {}, rightBCMatrix_?MatrixQ, {x_ /; ! NumericQ[x], xa_, xb_, xm_ : False},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
		Evans[{\[FormalLambda], \[FormalLambda]0}, AMatrix, N@SelectNegativeEigenvectors[AMatrix /. \[FormalLambda] -> \[FormalLambda]0, x], rightBCMatrix, {x, xa, xb, xm},opts]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, AMatrix_?MatrixQ, {}, {}, {x_ /; ! NumericQ[x], xa_, xb_, xm_ : False},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
		Evans[{\[FormalLambda], \[FormalLambda]0},
			AMatrix, N@SelectNegativeEigenvectors[AMatrix /. \[FormalLambda] -> \[FormalLambda]0, x], N@SelectPositiveEigenvectors[AMatrix /. \[FormalLambda] -> \[FormalLambda]0, x], {x, xa,xb, xm},opts]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, AMatrix_?MatrixQ, leftBCMatrix_?MatrixQ, {{}, {}}, {x_ /; ! NumericQ[x], xa_, xb_, xm_ : False},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
		Evans[{\[FormalLambda], \[FormalLambda]0},	AMatrix, leftBCMatrix, N@SelectPositiveEigenvectors[AMatrix /. \[FormalLambda] -> \[FormalLambda]0, x], {x, xa,	xb, xm},opts]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, {AMatrix_?MatrixQ, leftBCMatrix_?MatrixQ, rightBCMatrix_?MatrixQ, {x_ /; !NumericQ[x], xa_, xb_, xm_ : False}},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
    Evans[{\[FormalLambda], \[FormalLambda]0}, AMatrix, leftBCMatrix, rightBCMatrix, {x, xa, xb, xm},opts]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, AMatrix_?MatrixQ, leftBCMatrix_?MatrixQ,	rightBCMatrix_?MatrixQ, {x_ /; !NumericQ[x], xaa_, xbb_, xm_ : False},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
    Module[{len, subsets, newYs, leftYICs, rightYICs, phiLeftVector, phiRightVector, LeftBCSolution, RightBCSolution, yLeft, yRight,
	  phiLeft, phiRight, LeftPositiveEigenvalues, RightNegativeEigenvalues, phiLeftICs, phiRightICs, QQ, solutionFromRight,
	  solutionFromLeft, det, matchPoint,lenLeft,lenRight,subsetsLeft,subsetsRight,QLeft,QRight,xa,xb},

	If[!NumericQ[xaa], xa = (xaa /. \[FormalLambda] -> \[FormalLambda]0), xa = xaa];
	If[!NumericQ[xbb], xb = (xbb /. \[FormalLambda] -> \[FormalLambda]0), xb = xbb];

	If[!( xa < xb ),Message[Evans::incorrectRange,xa,xb];Return[$Failed]];

	If[Length[OptionValue[NormalizationConstants]]!=2,Message[Evans::normalizationConstants]; Return[$Failed]];

	If[NumericQ[xm /. \[FormalLambda] -> \[FormalLambda]0], matchPoint = xm /. \[FormalLambda] -> \[FormalLambda]0, matchPoint = (xb + xa) / 2];
  
  If[!(xa <= matchPoint <= xb), Message[Evans::matchPointIncorrect, xa, xb, matchPoint]; Return[$Failed]];

(* Check some conditions are true, square matrix with full rank and homogeneous BCs*)
  If[Length[AMatrix] != Length[Transpose[AMatrix]], Message[Evans::nonSquareMatrix];Return[$Failed]];

	(* Actually, don't necessarily need full rank, try without.
	If[Length[AMatrix] != MatrixRank[AMatrix], Message[Evans::matrixRank];Return[$Failed]]; *)

  (* Check that the eigenvalue appears in the matrix A, give warning otherwise*)
  If[MatrixQ[AMatrix, FreeQ[#, \[FormalLambda]] &],Message[Evans::noEigenvalue,AMatrix,\[FormalLambda]]];

(* Check that the matrix components are numerical, at least at the endpoints *)
  If[!MatrixQ[AMatrix /. \[FormalLambda] -> \[FormalLambda]0 /. x -> xa, NumericQ],Message[Evans::nonNumericalMatrix,AMatrix,xa];Return[$Failed]];
  If[!MatrixQ[AMatrix /. \[FormalLambda] -> \[FormalLambda]0 /. x -> xb, NumericQ],Message[Evans::nonNumericalMatrix,AMatrix,xb];Return[$Failed]];

	len = Length[AMatrix];
  newYs = Through[Array[\[FormalY], {len}][x]];

  (* Initial conditions for shooting from the LHS *)
  LeftBCSolution = Quiet@Solve[leftBCMatrix.newYs == 0, newYs];
	leftYICs = NullSpace[leftBCMatrix /. x -> xa /. \[FormalLambda] -> \[FormalLambda]0, Method -> "DivisionFreeRowReduction"];
  lenLeft = Length[leftYICs];
  subsetsLeft = Subsets[Range[len], {lenLeft}];

  (* Initial conditions for shooting from the RHS *)
  RightBCSolution = Quiet@Solve[rightBCMatrix.newYs == 0, newYs];
	rightYICs = NullSpace[rightBCMatrix /. x -> xb /. \[FormalLambda] -> \[FormalLambda]0, Method -> "DivisionFreeRowReduction"];
  lenRight = Length[rightYICs];
  subsetsRight = Subsets[Range[len], {lenRight}];

  (* Check the initial conditions for each side are enough *)
  If[Length[LeftBCSolution] != 1, Message[Evans::boundarySolutionFailed, xa];Return[$Failed]];
  If[Length[RightBCSolution] != 1, Message[Evans::boundarySolutionFailed, xb];Return[$Failed]];
  If[Length[leftYICs] + Length[rightYICs] != len, Message[Evans::boundaryConditionRank];Return[$Failed]];

  (* Generate two sets of Phi vaiables, these will be the matrix minors *)
  phiLeftVector = Table[\[FormalPhi]L[i][x], {i, 1, Length[subsetsLeft]}];
  phiRightVector = Table[\[FormalPhi]R[i][x], {i, 1, Length[subsetsRight]}];

  (* Full set of Initial Conditions for the left and right sides, with the BCs incorporated  *)
  yLeft = Transpose[leftYICs + Table[newYs, {lenLeft}] /. LeftBCSolution[[1]] /. Thread[newYs -> 0]];
  yRight = Transpose[rightYICs + Table[newYs, {lenRight}] /. RightBCSolution[[1]] /. Thread[newYs -> 0]];

  (* Use the initial conditions on the Y vectors to generate initial conditions for the minors phi *)

  phiLeft = (Det[(yLeft /. x -> xa /. \[FormalLambda] -> \[FormalLambda]0)[[#]]]& /@ subsetsLeft);
  phiRight = (Det[(yRight /. x -> xb /. \[FormalLambda] -> \[FormalLambda]0)[[#]]]& /@ subsetsRight);

  (* Find the exponentially growing modes from each side, positive or negative eigenvalues depending on which side *)

  LeftPositiveEigenvalues = Select[Eigenvalues[AMatrix /. x -> xa /. \[FormalLambda] -> \[FormalLambda]0], Re[#] > 0&];
  RightNegativeEigenvalues = Select[Eigenvalues[AMatrix /. x -> xb /. \[FormalLambda] -> \[FormalLambda]0], Re[#] < 0&];

  (* Apply the initial conditions for the left and right solutions *)

  phiLeftICs = Thread[Through[Array[\[FormalPhi]L, {Length[subsetsLeft]}][xa]] == phiLeft];
  phiRightICs = Thread[Through[Array[\[FormalPhi]R, {Length[subsetsRight]}][xb]] == phiRight];
   
  
  (* Calculate the Q matrix (phi' = Q phi) for each side *)
  QLeft = qMatrix[len, lenLeft] /. \[FormalCapitalA][i_, j_] :> AMatrix[[i, j]] /. \[FormalLambda] -> \[FormalLambda]0;
  QRight = qMatrix[len, lenRight] /. \[FormalCapitalA][i_, j_] :> AMatrix[[i, j]] /. \[FormalLambda] -> \[FormalLambda]0;


  (* Solve for integrating from the left and right *)
 
  solutionFromLeft = NDSolve[{Thread[D[phiLeftVector, x] == (QLeft -  	First[OptionValue[NormalizationConstants]] Total[Re@LeftPositiveEigenvalues] IdentityMatrix[Length[QLeft]]).phiLeftVector], phiLeftICs}, Array[\[FormalPhi]L, {Length[subsetsLeft]}], {x, xa, xb}, MaxStepFraction->0.01][[1]];
  solutionFromRight = NDSolve[{Thread[D[phiRightVector, x] == (QRight -  	Last[OptionValue[NormalizationConstants]] Total[Re@RightNegativeEigenvalues] IdentityMatrix[Length[QRight]]).phiRightVector], phiRightICs}, Array[\[FormalPhi]R, {Length[subsetsRight]}], {x, xa, xb},MaxStepFraction->0.01][[1]];

  (* Laplace Expanded Derivative of the determinant*)

  det = Total@Table[\[FormalPhi]L[i][x] \[FormalPhi]R[Complement[Range[len], i]][x] (-1)^(Total[Range[lenLeft] + i]) //. reprules /. Thread[subsetsLeft -> Range[Length[subsetsLeft]]] /. Thread[subsetsRight -> Range[Length[subsetsRight]]], {i, subsetsLeft}];
  
  (* Return the determinant, multiplied by the Evans factor to make it independent of the matching point *)
  
  Exp[-Integrate[Tr[AMatrix /. \[FormalLambda] -> \[FormalLambda]0], {x, xa, matchPoint}]] det /. x -> matchPoint /. solutionFromRight /. solutionFromLeft]

ToLinearMatrixForm[eqn_, bcs_?ListQ, depvars_, x_ /; (!NumericQ[x] && !ListQ[x])] :=
    If[bcs == {}, ToLinearMatrixForm[eqn, bcs, depvars, {x, 0, 0}], Message[ToLinearMatrixForm::noLimits];Return[$Failed]]

ToLinearMatrixForm[eqn_, BCs_?ListQ, depvars_, {x_ /; !NumericQ[x], xa_, xb_}] :=
    Module[{allVariables, nDepVars, allVariablesInd, originalYVariablesInd, nODEi, eqns, newYs, linearisedEqn, YDerivativeVector, FVector, originalYVariables, depVariables, epsilonEqn, newYSubs, 
	    newYDerivSubs, AMatrix, leftBCs, rightBCs, leftBCMatrix, rightBCMatrix, leftBCVector={}, rightBCVector={}, undifferentiatedVariables = {}, sol, undifsol},
	    
	    (* Cover for the case of a single variable or multiple: *)
	depVariables = Flatten[{depvars}];
  If[!ListQ[eqn], eqns = {eqn}, eqns = eqn];
	nDepVars = Length[depVariables];

		Table[
		(* Find highest derivative that occurs for each dependent variable   *)

			nODEi[i] = Max[0, Cases[eqns, Derivative[m_][depVariables[[i]]][x] :> m, Infinity]];
			
			(* Produce a list of all the derivatives for that variable *)
			allVariablesInd[i] = Table[D[depVariables[[i]][x], {x, n}], {n, 0, nODEi[i]}];

			(* If no derivatives exist, then that variable is undifferentiated and need to catch it. 
				 Else put all but the highest derivative into *)
			If[Length[allVariablesInd[i]] == 1,
				AppendTo[undifferentiatedVariables, allVariablesInd[i][[1]]];originalYVariablesInd[i] = {},
				originalYVariablesInd[i] = Most@allVariablesInd[i]],

			{i, nDepVars}];

		allVariables = Flatten[Table[allVariablesInd[i], {i, 1, nDepVars}]];
	
	(* Combine the solutions together *)
  originalYVariables = Flatten[Table[originalYVariablesInd[i], {i, 1, nDepVars}]];



  epsilonEqn = eqns /. Thread[allVariables -> (allVariables \[FormalE])];

	(* Linearise the equation in case it isn't already *)
	linearisedEqn = Normal@Series[epsilonEqn, {\[FormalE], 0, 1}];

	(* Check if the linearised equation is the same as the original, if not throw a warning to make it explicit *)
	If[! (Thread[(((linearisedEqn /. Equal -> Subtract) - (epsilonEqn /.Equal -> Subtract) // Simplify) == Table[0, nDepVars])] ===	True),
		Message[ToLinearMatrixForm::linearised, linearisedEqn /. \[FormalE] -> 1]];

	(* Replace all the original variables with a set indexed by Y *)
  newYs = Through[Array[\[FormalY], {Length[originalYVariables]}][x]];
  newYSubs = Thread[originalYVariables -> newYs];
  newYDerivSubs = D[newYSubs, x];

	(* Solve for the derivatives of Y *)
	
  sol = Solve[Flatten@Join[Select[Thread[D[originalYVariables /. newYSubs, x] == D[originalYVariables, x] /. newYSubs], 
	        FreeQ[#, Alternatives @@ depVariables]&], {linearisedEqn /. \[FormalE] -> 1} //. newYSubs /. newYDerivSubs], 
	        Join[D[originalYVariables /. newYSubs, x], undifferentiatedVariables]];
	
	(* Check that a single solution has been found *)
  If[Length[sol] != 1, Message[ToLinearMatrixForm::somethingWrong];Return[$Failed]];
 
	
	(* Find the right hand side of Y' = A Y + F*)
	YDerivativeVector =	D[originalYVariables /. newYSubs, x] /.sol[[1]];
	AMatrix = Coefficient[#, newYs]& /@ YDerivativeVector;
	FVector = YDerivativeVector /. Thread[newYs -> 0];

	(* If the equation is inhomogeneous, print warning to say the CMM doesn't work *)
  If[AnyTrue[FVector, (!PossibleZeroQ[#])&], Message[ToLinearMatrixForm::inhomogeneous]];

	(*
	Don't think this is doing anything...
	If[Length[undifferentiatedVariables] > 0, undifsol = (DSolve[Take[sol[[1]], -Length[undifferentiatedVariables]] /. Rule -> Equal, depvars, x] // Quiet), undifsol = {}];
	
	*)
	(*Generate the boundary condition matrices and vectors *)
	leftBCs = Select[BCs /. (Thread[originalYVariables -> newYs] /. x -> xa) /. (sol[[1]] /. x -> xa), !FreeQ[#, \[FormalY]]&];
  rightBCs = Select[BCs /. (Thread[originalYVariables -> newYs] /. x -> xb) /. (sol[[1]] /. x -> xb), !FreeQ[#, \[FormalY]]&];

	If[Length[leftBCs] > 0,
		leftBCMatrix = Coefficient[#, newYs /. x -> xa] & /@ (leftBCs /. Equal -> Subtract);
		leftBCVector = -(leftBCs /. Equal -> Subtract) /.	Thread[(newYs /. x -> xa) -> 0],
		  leftBCMatrix = {};
	];

	If[Length[rightBCs] > 0,
		rightBCMatrix =	Coefficient[#, newYs /. x -> xb] & /@ (rightBCs /.	Equal -> Subtract);
		rightBCVector = -(rightBCs /. Equal -> Subtract) /. Thread[(newYs /. x -> xb) -> 0];,
		   rightBCMatrix = {};
	];
	
	(* Return the solutions *)
	If[AnyTrue[Join[leftBCVector,rightBCVector],(!PossibleZeroQ[#])&],
		Message[ToLinearMatrixForm::inhomogeneousBCs,leftBCVector,rightBCVector];
		{AMatrix, {leftBCMatrix, leftBCVector}, {rightBCMatrix, rightBCVector}, {x, xa, xb}},

	    {AMatrix, leftBCMatrix, rightBCMatrix, {x, xa, xb}}
	]
]

ToLinearMatrixForm[eqns_?ListQ, BCs_?ListQ, {depvarLeft_, depvarRight_}, {x_, xa_, xmatch_, xb_}] :=

(* Currently only works if the two dependent variables are given in that specific order *)

		Module[{leftEqns, rightEqns, leftBCs, rightBCs, FMatrix, GMatrix,	leftAMatrix, rightAMatrix, leftBCMatrix, rightBCMatrix, stuff, xLeft, xRight,interfaceBCs,flatBCs},
			flatBCs = Flatten[BCs];
			leftEqns = Select[eqns, ! FreeQ[#, depvarLeft] &];
			rightEqns = Select[eqns, ! FreeQ[#, depvarRight] &];
			leftBCs =		Select[flatBCs, !	FreeQ[#, \[FormalA]_[xa] /; ! (\[FormalA] === Derivative),	All] &];
			rightBCs =	Select[flatBCs, !	FreeQ[#, \[FormalA]_[xb] /; ! (\[FormalA] === Derivative),	All] &];
			interfaceBCs =	Select[flatBCs,   !	FreeQ[#, \[FormalA]_[xmatch] /; ! (\[FormalA] === Derivative),	All] &];
			FMatrix = ExtractInterface[interfaceBCs, depvarLeft, {x, xmatch}];
			GMatrix = ExtractInterface[interfaceBCs, depvarRight, {x, xmatch}];
			{leftAMatrix, leftBCMatrix, stuff, xLeft} =		  ToLinearMatrixForm[leftEqns, leftBCs, depvarLeft, {x, xa, xmatch}];
			{rightAMatrix, stuff, rightBCMatrix, xRight} = 	ToLinearMatrixForm[rightEqns, rightBCs,	depvarRight, {x, xmatch, xb}];
			(*Return the system for using in CMM function*)
			{{leftAMatrix, rightAMatrix}, leftBCMatrix,rightBCMatrix, {FMatrix,	GMatrix}, {x, xa, xmatch, xb}}

		]

ExtractInterface[eqn_, depvars_, {x_,xmatch_}] :=
		Module[{undifferentiatedVariables = {},neweqn,depVariables,nDepVars,nODEi,eqns,allVariablesInd,originalYVariablesInd,allVariables,originalYVariables,newYs,newYSubs},
		(*Cover for the case of a single variable or multiple:*)
			neweqn = eqn /. \[FormalY]_[xmatch] -> \[FormalY][x];
			depVariables = Flatten[{depvars}];
			If[! ListQ[neweqn], eqns = {neweqn}, eqns = neweqn];
			nDepVars = Length[depVariables];
			Table[(*Find highest derivative that occurs for each dependent variable*)
				nODEi[i] =	Max[0, Cases[eqns, Derivative[m_][depVariables[[i]]][a_] :> m, Infinity]];
				(*Produce a list of all the derivatives for that variable*)
				allVariablesInd[i] =	Table[D[depVariables[[i]][x], {x, n}], {n, 0, nODEi[i]}];
				(*If no derivatives exist,	then that variable is undifferentiated and need to catch it. Else put all in *)
				If[Length[allVariablesInd[i]] == 1,
					AppendTo[undifferentiatedVariables, allVariablesInd[i][[1]]];
					originalYVariablesInd[i] = {},
					originalYVariablesInd[i] = allVariablesInd[i]], {i, nDepVars}];
			allVariables = Flatten[Table[allVariablesInd[i], {i, 1, nDepVars}]];
			(*Combine the solutions together*)
			originalYVariables =	Flatten[Table[originalYVariablesInd[i], {i, 1, nDepVars}]];

			(*Replace all the original variables with a set indexed by Y*)
			newYs = Through[Array[\[FormalY], {Length[originalYVariables]}][x]];
			newYSubs = Thread[originalYVariables -> newYs];
			Transpose[Table[Coefficient[neweqn /. newYSubs /. Equal -> Subtract, i], {i, newYs}]]
		]

SelectNegativeEigenvectors[mat_?MatrixQ, x_ /; !NumericQ[x]] := Module[{limitMatrix,eigvecs},

		limitMatrix = Limit[mat, x -> -\[Infinity]];
		If[!MatrixQ[limitMatrix, NumericQ], Message[SelectNegativeEigenvectors::nonNumericalMatrix,limitMatrix];Return[$Failed]];
		eigvecs=Extract[#[[2]], Position[#[[1]], a_ /; Re[a] < 0, 1]] &@	Eigensystem[limitMatrix];
	  If[Length[eigvecs]==0, Message[SelectNegativeEigenvectors::noNegativeEigenvectors];Return[$Failed],Return[eigvecs]]

	]

SelectPositiveEigenvectors[mat_?MatrixQ, x_ /; !NumericQ[x]] := Module[{limitMatrix,eigvecs},
		limitMatrix=Limit[mat, x -> \[Infinity]];
		If[!MatrixQ[limitMatrix, NumericQ], Message[SelectPositiveEigenvectors::nonNumericalMatrix,limitMatrix];Return[$Failed]];
		eigvecs=Extract[#[[2]], Position[#[[1]], a_ /; Re[a] > 0, 1]] &@	Eigensystem[limitMatrix];
	  If[Length[eigvecs]==0, Message[SelectPositiveEigenvectors::noPositiveEigenvectors];Return[$Failed],Return[eigvecs]]

	]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, {{ALeftMatrix_?MatrixQ,	ARightMatrix_?MatrixQ},
	leftBCMatrix_?MatrixQ, rightBCMatrix_?MatrixQ, {FMatrix_?MatrixQ,	GMatrix_?MatrixQ}, {x_ /; ! NumericQ[x], xa_, xm_, xb_}},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
		Evans[{\[FormalLambda], \[FormalLambda]0}, {ALeftMatrix,	ARightMatrix}, leftBCMatrix, rightBCMatrix,	{FMatrix,	GMatrix}, {x, xa, xm, xb},opts]

Evans[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?
		NumericQ}, {ALeftMatrix_?MatrixQ,	ARightMatrix_?MatrixQ}, leftBCMatrix_?MatrixQ, rightBCMatrix_?MatrixQ,
	{FMatrix_?MatrixQ,	GMatrix_?MatrixQ}, {x_ /; !NumericQ[x], xaa_, xmm_, xbb_},opts:OptionsPattern[{NormalizationConstants -> {1, 1}}]] :=
		Module[{dettt, len, subsets, newYs, leftYICs, rightYICs, phiLeftVector, phiRightVector, LeftBCSolution, RightBCSolution, yLeft, yRight,
			phiLeft, phiRight, LeftPositiveEigenvalues, RightNegativeEigenvalues, phiLeftICs, phiRightICs, QQ, solutionFromRight,
			solutionFromLeft, det, matchPoint,lenLeft,lenRight,subsetsLeft,subsetsRight,QLeft,QRight,xa,xb,xm},

			If[! NumericQ[xaa], xa = (xaa /. \[FormalLambda] -> \[FormalLambda]0),	xa = xaa];
			If[! NumericQ[xmm], xm = (xmm /. \[FormalLambda] -> \[FormalLambda]0),	xm = xmm];
			If[! NumericQ[xbb],	xb = (xbb /. \[FormalLambda] -> \[FormalLambda]0),  xb = xbb];

			If[Length[OptionValue[NormalizationConstants]]!=2,Message[Evans::normalizationConstants]; Return[$Failed]];

		(*Check some conditions are true *)
		If[Length[ALeftMatrix] != Length[Transpose[ALeftMatrix]],
			Message[Evans::nonSquareMatrix]; Return[$Failed]];
		If[Length[ARightMatrix] != Length[Transpose[ARightMatrix]],
			Message[Evans::nonSquareMatrix]; Return[$Failed]];
		(*Check that the eigenvalue appears in the matrix A*)
		If[MatrixQ[ALeftMatrix, FreeQ[#, \[FormalLambda]] &] &&
				MatrixQ[ARightMatrix, FreeQ[#, \[FormalLambda]] &],
			Message[Evans::noEigenvalue,
				{ALeftMatrix,ARightMatrix}, \[FormalLambda]];];
		(*Check that the matrix components are numerical,	at least at the endpoints*)
		If[! MatrixQ[ALeftMatrix /. x -> xa /. \[FormalLambda] -> \[FormalLambda]0 ,NumericQ],
			Message[Evans::nonNumericalMatrix, ALeftMatrix, xa];	Return[$Failed]];
		If[! MatrixQ[ARightMatrix /. x -> xa /. \[FormalLambda] -> \[FormalLambda]0 ,NumericQ],
			Message[Evans::nonNumericalMatrix, ARightMatrix, xa]; Return[$Failed]];
		If[! MatrixQ[ALeftMatrix /. x -> xb /. \[FormalLambda] -> \[FormalLambda]0 ,NumericQ],
			Message[Evans::nonNumericalMatrix, ALeftMatrix, xb];	Return[$Failed]];
		If[! MatrixQ[ARightMatrix /. x -> xb /. \[FormalLambda] -> \[FormalLambda]0 ,NumericQ],
			Message[Evans::nonNumericalMatrix, ARightMatrix,	xb]; Return[$Failed]];


		If[Length[ARightMatrix] != Length[ALeftMatrix],
			Message[Evans::MatrixSizesDiffer, ALeftMatrix, ARightMatrix];Return[$Failed]];

		len = Length[ARightMatrix];

		If[len>10,	Message[Evans::InterfaceTooBig, len];Return[$Failed]];


		newYs = Through[Array[\[FormalY], {len}][x]];

		(*Initial conditions for shooting from the LHS*)
		LeftBCSolution =
				Quiet@Solve[leftBCMatrix.newYs == 0, newYs];
		leftYICs =
				NullSpace[leftBCMatrix /. x -> xa /. \[FormalLambda] -> \[FormalLambda]0,	Method -> "DivisionFreeRowReduction"];
		lenLeft = Length[leftYICs];
		subsetsLeft = Subsets[Range[len], {lenLeft}];
		(*Initial conditions for shooting from the RHS*)
		RightBCSolution =
				Quiet@Solve[rightBCMatrix.newYs == 0, newYs];
		rightYICs =
				NullSpace[rightBCMatrix /. x -> xb /. \[FormalLambda] -> \[FormalLambda]0,	Method -> "DivisionFreeRowReduction"];
		lenRight = Length[rightYICs];
		subsetsRight = Subsets[Range[len], {lenRight}];
		(*Check the initial conditions for each side are enough*)
		If[Length[LeftBCSolution] != 1,
			Message[Evans::boundarySolutionFailed, xa];
			Return[$Failed]];
		If[Length[RightBCSolution] != 1,
			Message[Evans::boundarySolutionFailed, xb];
			Return[$Failed]];
		If[Length[leftYICs] + Length[rightYICs] != len,
			Message[Evans::boundaryConditionRank];
			Return[$Failed]];
		(*Generate two sets of Phi vaiables,these will be the matrix minors*)

		phiLeftVector =
				Table[\[FormalPhi]L[i][x], {i, 1, Length[subsetsLeft]}];
		phiRightVector =
				Table[\[FormalPhi]R[i][x], {i, 1, Length[subsetsRight]}];
		(*Full set of Initial Conditions for the left and right sides,	with the BCs incorporated*)
		yLeft =
				Transpose[leftYICs + Table[newYs, {lenLeft}] /. LeftBCSolution[[1]] /.Thread[newYs -> 0]];
		yRight =
				Transpose[rightYICs + Table[newYs, {lenRight}] /. RightBCSolution[[1]] /.Thread[newYs -> 0]];
		(*Use the initial conditions on the Y vectors to generate initial conditions for the minors phi*)
		phiLeft = (Det[(yLeft /.   x -> xa /. \[FormalLambda] -> \[FormalLambda]0)[[#]]] & /@subsetsLeft);
		phiRight = (Det[(yRight /. x -> xb /. \[FormalLambda] -> \[FormalLambda]0)[[#]]] & /@subsetsRight);
		(*Find the exponentially growing modes from each side,	positive or negative eigenvalues depending on which side*)
		LeftPositiveEigenvalues =
				Select[Eigenvalues[
					ALeftMatrix /. x -> xa /. \[FormalLambda] -> \[FormalLambda]0],
					Re[#] > 0 &];
		RightNegativeEigenvalues =
				Select[Eigenvalues[
					ARightMatrix /. x -> xb /. \[FormalLambda] -> \[FormalLambda]0],
					Re[#] < 0 &];
		(*Apply the initial conditions for the left and right solutions*)
		phiLeftICs =
				Thread[Through[Array[\[FormalPhi]L, {Length[subsetsLeft]}][xa]] ==
						phiLeft];
		phiRightICs =
				Thread[Through[Array[\[FormalPhi]R, {Length[subsetsRight]}][xb]] ==
						phiRight];
		(*Calculate the Q matrix (phi'=Q phi) for each side*)
		QLeft =
				qMatrix[len, lenLeft] /. \[FormalCapitalA][i_, j_] :>
						ALeftMatrix[[i, j]] /. \[FormalLambda] -> \[FormalLambda]0;
		QRight =
				qMatrix[len, lenRight] /. \[FormalCapitalA][i_, j_] :>
						ARightMatrix[[i, j]] /. \[FormalLambda] -> \[FormalLambda]0;

		(*Solve for integrating from the left and right*)
		solutionFromLeft =
				NDSolve[{Thread[
					D[phiLeftVector, x] == (QLeft -	First[OptionValue[NormalizationConstants]] Total[Re@LeftPositiveEigenvalues] IdentityMatrix[Length[QLeft]]).phiLeftVector], phiLeftICs},
					Array[\[FormalPhi]L, {Length[subsetsLeft]}], {x, xa, xm},
					MaxStepFraction -> 0.05][[1]];
		solutionFromRight =
				NDSolve[{Thread[
					D[phiRightVector,	x] == (QRight -	Last[OptionValue[NormalizationConstants]] Total[Re@RightNegativeEigenvalues] IdentityMatrix[Length[QRight]]).phiRightVector], phiRightICs},
					Array[\[FormalPhi]R, {Length[subsetsRight]}], {x, xm, xb},
					MaxStepFraction -> 0.05][[1]];

		(* Now we need to account for the jump conditions at the interface, so instead of the normal determinant it needs
		   modifying by multiplication by the matrices F and G. *)

		det = Total@Table[\[FormalPhi]\[FormalPhi]L[i][x]
					\[FormalPhi]\[FormalPhi]R[Complement[Range[len], i]][x] (-1)^(Total[Range[lenLeft] + i]), {i, subsetsLeft}];

		dettt =
				det /.rulesFG[len]/.reprules2 /.Thread[subsetsLeft -> Range[Length[subsetsLeft]]]
						/.Thread[subsetsRight -> Range[Length[subsetsRight]]];

		Exp[-Integrate[
			Tr[ALeftMatrix /. \[FormalLambda] -> \[FormalLambda]0], {x, xa, xm}]]
			dettt /. {F[i_, j_] :> FMatrix[[i, j]],	G[i_, j_] :> GMatrix[[i, j]]} /.
				x -> xm /. \[FormalLambda] -> \[FormalLambda]0 /.	solutionFromRight /. solutionFromLeft
		]



End[]; (* End Private Context *)

EndPackage[];
