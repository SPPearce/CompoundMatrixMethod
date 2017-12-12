(* ::Package:: *)

(* Mathematica Package  *)

(* :Title: CompoundMatrixMethod *)
(* :Author: Simon Pearce <simon.pearce@manchester.ac.uk> *)
(* :Context: CompoundMatrixMethod` *)
(* :Version: 0.3 *)
(* :Date: 2017-12-07 *)

(* :Mathematica Version: 9+ *)
(* :Copyright: (c) 2017 Simon Pearce *)

BeginPackage["CompoundMatrixMethod`"];

Unprotect["CompoundMatrixMethod`*"];

CompoundMatrixMethod::usage = "\
CompoundMatrixMethod[{k,k0},sys] evaluates the Evans function from the Compound Matrix Method with potential eigenvalue k=k0, for the system defined from ToLinearMatrixForm.
CompoundMatrixMethod[{k,k0},A,B,C,{x,x0,x1}] evaluates the Evans function from the Compound Matrix Method with k=k0. Here the linear matrix ODE is given by dy/dx=A.y, with B.y=0 at x=x0 and C.y=0 at x=x1.
For either case a complex number is returned, zeroes of this function correspond to zeroes of the original eigenvalue equation.";

ToLinearMatrixForm::usage = "\
ToLinearMatrixForm[eqn,{},depvars,x] takes a list of differential equations in the dependent variables depvars and independent variable x and puts the equations into linear matrix form dy/dx = A.y
ToLinearMatrixForm[eqn,bcs,depvars,{x,x0,x1}] also includes the boundary conditions evaluated at x=x0 and x=x1.";

CompoundMatrixMethod::boundarySolutionFailed = "Applying the boundary conditions at `1` failed, perhaps you don't have any conditions there.";
CompoundMatrixMethod::nonSquareMatrix = "The matrix A is not square.";
CompoundMatrixMethod::matrixRank = "The matrix A has rank less than its dimensions.";

CompoundMatrixMethod::boundaryConditionRank = "The rank of the boundary condition matrices does not sum to the size of A";
CompoundMatrixMethod::nonZeroBoundaryConditions = "Non-zero boundary condition given, "

CompoundMatrixMethod::nonNumericalMatrix = "Matrix `1` is not numerical at x = `2` "
CompoundMatrixMethod::noEigenvalue = "Matrix `1` does not contain the potential eigenvalue `2` "


ToLinearMatrixForm::somethingWrong = "Something has gone wrong!";

ToLinearMatrixForm::matrixOnly = "Incorrect number of boundary conditions given (`1` compared to matrix dimension `2`), returning the matrix only";
ToLinearMatrixForm::inhomogeneous = "Inhomogeneous equation, does not work with the Compound Matrix Method";
ToLinearMatrixForm::linearised = "Original Equation not linear, linearized equation given by `1`";
ToLinearMatrixForm::noLimits = "Please supply limits for the independent variable when providing boundary conditions";

SelectPositiveEigenvectors::usage = "\
SelectPositiveEigenvectors[matrix, x] selects the eigenvectors which correspond to negative real part as x->-inf"

SelectNegativeEigenvectors::usage = "\
SelectNegativeEigenvectors[matrix] selects the eigenvectors which correspond to positive real part at  x->inf"

SelectPositiveEigenvectors::nonNumericalMatrix = "Matrix `1` is not numerical"
SelectNegativeEigenvectors::nonNumericalMatrix = "Matrix `1` is not numerical"

Begin["`Private`"]; (* Begin Private Context *)

(* Replacement rule to sort the lists of indices *)
reprules = \[FormalPhi][a_List] :> Signature[a] \[FormalPhi][Sort[a]];

(* Generation of the derivatives of the matrix minors, looping over  rule to sort the lists of indices *)
minorsDerivs[list_?VectorQ, len_?NumericQ] := minorsDerivs[list, len] =
    Sum[
	      Sum[\[FormalCapitalA][y, z] \[FormalPhi][list /. y -> z], {z, Union[Complement[Range[len], list], {y}]}],
	    {y, list}] /. reprules

qMatrix[len_?NumericQ, len2_?NumericQ] := qMatrix[len, len2] =
    Module[{minorsTab},
  minorsTab = Table[minorsDerivs[ii, len], {ii, Subsets[Range[len], {len2}]}]
		  /. Thread[Subsets[Range[len], {len2}] -> Range[Length[Subsets[Range[len], {len2}]]]];
  Transpose[Table[Coefficient[minorsTab, \[FormalPhi][i]], {i, 1, Binomial[len, len2]}]]
]
CompoundMatrixMethod[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, {AMatrix_?MatrixQ, leftBCMatrix_?MatrixQ, rightBCMatrix_?MatrixQ, {x_ /; !NumericQ[x], xa_?NumericQ, xb_?NumericQ, xm_ : False}}] :=
    CompoundMatrixMethod[{\[FormalLambda], \[FormalLambda]0}, AMatrix, {leftBCMatrix, 0}, {rightBCMatrix, 0}, {x, xa, xb, xm}]
CompoundMatrixMethod[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, {AMatrix_?MatrixQ, {leftBCMatrix_?MatrixQ, leftBCVector_}, {rightBCMatrix_?MatrixQ, rightBCVector_}, {x_ /; !NumericQ[x], xa_?NumericQ, xb_?NumericQ, xm_ : False}}] :=
    CompoundMatrixMethod[{\[FormalLambda], \[FormalLambda]0}, AMatrix, {leftBCMatrix, leftBCVector}, {rightBCMatrix, rightBCVector}, {x, xa, xb, xm}]
CompoundMatrixMethod[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, AMatrix_?MatrixQ, leftBCMatrix_?MatrixQ, rightBCMatrix_?MatrixQ, {x_ /; !NumericQ[x], xa_?NumericQ, xb_?NumericQ, xm_ : False}] :=
    CompoundMatrixMethod[{\[FormalLambda], \[FormalLambda]0}, AMatrix, {leftBCMatrix, 0}, {rightBCMatrix, 0}, {x, xa, xb, xm}]

CompoundMatrixMethod[{\[FormalLambda]_ /; !NumericQ[\[FormalLambda]], \[FormalLambda]0_?NumericQ}, AMatrix_?MatrixQ, {leftBCMatrix_?MatrixQ, leftBCVector_}, {rightBCMatrix_?MatrixQ, rightBCVector_}, {x_ /; !NumericQ[x], xa_?NumericQ, xb_?NumericQ, xm_ : False}] := Module[
  {len, subsets, newYs, leftYICs, rightYICs, phiLeftVector, phiRightVector, LeftBCSolution, RightBCSolution, yLeft, yRight,
	  phiLeft, phiRight, LeftPositiveEigenvalues, RightNegativeEigenvalues, phiLeftICs, phiRightICs, QQ, solutionFromRight,
	  solutionFromLeft, det, matchPoint,lenLeft,lenRight,subsetsLeft,subsetsRight,QLeft,QRight},
  If[(xa <= xm <= xb && NumericQ[xm]), matchPoint = xm, matchPoint = (xb + xa) / 2];

  (* Check some conditions are true, square matrix with full rank and homogeneous BCs*)
  If[Length[AMatrix] != Length[Transpose[AMatrix]], Message[CompoundMatrixMethod::nonSquareMatrix];Return[$Failed]];
  If[Length[AMatrix] != MatrixRank[AMatrix], Message[CompoundMatrixMethod::matrixRank];Return[$Failed]];
  (* Check the equations are inhomogeneous *)
  If[Max@Abs@leftBCVector != 0 || Max@Abs@rightBCVector != 0, Message[CompoundMatrixMethod::nonZeroBoundaryConditions];Return[$Failed]];

  (* Check that the eigenvalue appears in the matrix A *)
  If[MatrixQ[AMatrix, FreeQ[#, \[FormalLambda]] &],Message[CompoundMatrixMethod::noEigenvalue,AMatrix,\[FormalLambda]];Return[$Failed]];

  (* Check that the matrix components are numerical, at least at the endpoints *)
  If[!MatrixQ[AMatrix /. \[FormalLambda] -> \[FormalLambda]0 /. x -> xa, NumericQ],Message[CompoundMatrixMethod::nonNumericalMatrix,AMatrix,xa];Return[$Failed]];
  If[!MatrixQ[AMatrix /. \[FormalLambda] -> \[FormalLambda]0 /. x -> xb, NumericQ],Message[CompoundMatrixMethod::nonNumericalMatrix,AMatrix,xb];Return[$Failed]];

	len = Length[AMatrix];
  newYs = Through[Array[\[FormalY], {len}][x]];

  (* Initial conditions for shooting from the LHS *)
  LeftBCSolution = Quiet@Solve[leftBCMatrix.newYs == leftBCVector, newYs];
	leftYICs = NullSpace[leftBCMatrix /. x -> xa];
  lenLeft = Length[leftYICs];
  subsetsLeft = Subsets[Range[len], {lenLeft}];


  (* Initial conditions for shooting from the RHS *)
  RightBCSolution = Quiet@Solve[rightBCMatrix.newYs == rightBCVector, newYs];
	rightYICs = NullSpace[rightBCMatrix /. x -> xb];
  lenRight = Length[rightYICs];
  subsetsRight = Subsets[Range[len], {lenRight}];

  (* Check the initial conditions for each side are enough *)
  If[Length[LeftBCSolution] != 1, Message[CompoundMatrixMethod::boundarySolutionFailed, xa];Return[$Failed]];
  If[Length[RightBCSolution] != 1, Message[CompoundMatrixMethod::boundarySolutionFailed, xb];Return[$Failed]];
  If[Length[leftYICs] + Length[rightYICs] != len, Message[CompoundMatrixMethod::boundaryConditionRank];Return[$Failed]];

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
 
  solutionFromLeft = NDSolve[{Thread[D[phiLeftVector, x] == (QLeft - Total[Re@LeftPositiveEigenvalues] IdentityMatrix[Length[QLeft]]).phiLeftVector], phiLeftICs}, Array[\[FormalPhi]L, {Length[subsetsLeft]}], {x, xa, xb}][[1]];
  solutionFromRight = NDSolve[{Thread[D[phiRightVector, x] == (QRight - Total[Re@RightNegativeEigenvalues] IdentityMatrix[Length[QRight]]).phiRightVector], phiRightICs}, Array[\[FormalPhi]R, {Length[subsetsRight]}], {x, xa, xb}][[1]];

  (* Laplace Expanded Derivative of the determinant*)

  det = Total@Table[\[FormalPhi]L[i][x] \[FormalPhi]R[Complement[Range[len], i]][x] (-1)^(Total[Range[lenLeft] + i]) //. reprules /. Thread[subsetsLeft -> Range[Length[subsetsLeft]]] /. Thread[subsetsRight -> Range[Length[subsetsRight]]], {i, subsetsLeft}];
  
  (* Return the determinant, multiplied by the Evans factor to make it independent of the matching point *)
  
  Exp[-Integrate[Tr[AMatrix /. \[FormalLambda] -> \[FormalLambda]0], {x, xa, matchPoint}]] det /. x -> matchPoint /. solutionFromRight /. solutionFromLeft]

ToLinearMatrixForm[eqn_, bcs_?ListQ, depvars_, x_] := 
    If[bcs == {}, ToLinearMatrixForm[eqn, bcs, depvars, {x, 0, 0}], Message[ToLinearMatrixForm::noLimits];Return[$Failed]]

ToLinearMatrixForm[eqn_, BCs_?ListQ, depvars_, {x_, xa_, xb_}] := 
    Module[{allVariables, nDepVars, allVariablesInd, originalYVariablesInd, nODEi, eqns, newYs, linearisedEqn, YDerivativeVector, FVector, originalYVariables, depVariables, epsilonEqn, newYSubs, 
	    newYDerivSubs, AMatrix, leftBCs, rightBCs, leftBCMatrix, rightBCMatrix, leftBCVector, rightBCVector, undifferentiatedVariables = {}, sol, undifsol},
	    
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

	(* If the equation is inhomogeneous, print message to say the CMM doesn't work *)
  If[AnyTrue[FVector, (!PossibleZeroQ[#])&], Message[ToLinearMatrixForm::inhomogeneous]];

	(* If no BCs were given, just return the Matrix and give a warning *)
  If[Length[BCs] != Length[AMatrix], Message[ToLinearMatrixForm::matrixOnly,Length[BCs],Length[AMatrix]];Return[AMatrix]];

	(*
	Don't think this is doing anything...
	If[Length[undifferentiatedVariables] > 0, undifsol = (DSolve[Take[sol[[1]], -Length[undifferentiatedVariables]] /. Rule -> Equal, depvars, x] // Quiet), undifsol = {}];
	
	*)
	(*Generate the boundary condition matrices and vectors *)
	leftBCs = Select[BCs /. (Thread[originalYVariables -> newYs] /. x -> xa) /. (sol[[1]] /. x -> xa), !FreeQ[#, \[FormalY]]&];
  rightBCs = Select[BCs /. (Thread[originalYVariables -> newYs] /. x -> xb) /. (sol[[1]] /. x -> xb), !FreeQ[#, \[FormalY]]&];
  leftBCMatrix = Coefficient[#, newYs /. x -> xa]& /@ (leftBCs /. Equal -> Subtract);
  rightBCMatrix = Coefficient[#, newYs /. x -> xb]& /@ (rightBCs /. Equal -> Subtract);
  leftBCVector = -(leftBCs /. Equal -> Subtract) /. Thread[(newYs /. x -> xa) -> 0];
  rightBCVector = -(rightBCs /. Equal -> Subtract) /. Thread[(newYs /. x -> xb) -> 0];
	
	(* Return the solutions *)
  {AMatrix, {leftBCMatrix, leftBCVector}, {rightBCMatrix, rightBCVector}, {x, xa, xb}}
]

SelectNegativeEigenvectors[mat_?MatrixQ, x_ /; !NumericQ[x]] := Module[{limitMatrix},

		limitMatrix = Limit[mat, x -> -\[Infinity]];
		If[!MatrixQ[limitMatrix, NumericQ], Message[SelectNegativeEigenvectors::nonNumericalMatrix,limitMatrix];Return[$Failed]];
		Extract[#[[2]], Position[#[[1]], a_ /; Re[a] < 0, 1]] &@	Eigensystem[limitMatrix]

	]

SelectPositiveEigenvectors[mat_?MatrixQ, x_ /; !NumericQ[x]] := Module[{limitMatrix},
		limitMatrix=Limit[mat, x -> \[Infinity]];
		If[!MatrixQ[limitMatrix, NumericQ], Message[SelectNegativeEigenvectors::nonNumericalMatrix,limitMatrix];Return[$Failed]];
		Extract[#[[2]], Position[#[[1]], a_ /; Re[a] > 0, 1]] &@	Eigensystem[limitMatrix]
	]

End[]; (* End Private Context *)

EndPackage[];
