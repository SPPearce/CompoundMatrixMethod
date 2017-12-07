(* ::Package:: *)

(* Mathematica Package  *)

(* :Title: CompoundMatrixMethod *)
(* :Author: Simon Pearce <simon.pearce@manchester.ac.uk> *)
(* :Context: CompoundMatrixMethod` *)
(* :Version: 0.1 *)
(* :Date: 2017-12-07 *)

(* :Mathematica Version: 9+ *)
(* :Copyright: (c) 2017 Simon Pearce *)

BeginPackage["CompoundMatrixMethod`"];

Unprotect["CompoundMatrixMethod`*"];

CompoundMatrixMethod::usage = "\
CompoundMatrixMethod[{k,kval},sys] evaluates the Evans function from the Compound Matrix Method with potential eigenvalue k=kval, for the system defined from ToLinearMatrixForm.
CompoundMatrixMethod[{k,kval},A,B,C,{x,x0,x1}] evaluates the Evans function from the Compound Matrix Method with k=kval. Here the linear matrix ODE is given by dy/dx=A.y, with B.y=0 at x=x0 and C.y=0 at x=x1.
For either case a complex number is returned, zeroes of this function correspond to zeroes of the original eigenvalue equation.";

ToLinearMatrixForm::usage = "\
ToLinearMatrixForm[eqn,{},depvars,x] takes a list of differential equations in the dependent variables depvars and independent variable x and puts the equations into linear matrix form dy/dx = A.y
ToLinearMatrixForm[eqn,bcs,depvars,{x,x0,x1}] also includes the boundary conditions evaluated at x=x0 and x=x1.";

Begin["`Private`"]; (* Begin Private Context *)

reprules=\[FormalPhi][a_List]:>Signature[a] \[FormalPhi][Sort[a]];
minorsDerivs[list_?VectorQ,len_]:=minorsDerivs[list,len]=Sum[Sum[\[FormalCapitalA][y,z] \[FormalPhi][list/.y->z],{z,Union[Complement[Range[len],list],{y}]}],{y,list}]/.reprules
qMatrix[len_?EvenQ,len2_]:=qMatrix[len,len2]=Module[{minorsTab},
minorsTab=Table[minorsDerivs[ii,len],{ii,Subsets[Range[len],{len2}]}]/.Thread[Subsets[Range[len],{len2}]->Range[Length[Subsets[Range[len],{len2}]]]];
Transpose[Table[Coefficient[minorsTab,\[FormalPhi][i]],{i,1,Binomial[len,len2]}]]
]
CompoundMatrixMethod[{\[Lambda]_/;!NumericQ[\[Lambda]],\[Lambda]\[Lambda]_?NumericQ},{Amat_?MatrixQ,bmat_,cmat_,{x_/;!NumericQ[x],xa_?NumericQ,xb_?NumericQ,xmatch_:False}},dampfactor_:1]:=CompoundMatrixMethod[{\[Lambda],\[Lambda]\[Lambda]},Amat,{bmat,0},{cmat,0},{x,xa,xb,xmatch}]
CompoundMatrixMethod[{\[Lambda]_/;!NumericQ[\[Lambda]],\[Lambda]\[Lambda]_?NumericQ},{Amat_?MatrixQ,{bmat_?MatrixQ,bvec_},{cmat_?MatrixQ,cvec_},{x_/;!NumericQ[x],xa_?NumericQ,xb_?NumericQ,xmatch_:False}},dampfactor_:1]:=CompoundMatrixMethod[{\[Lambda],\[Lambda]\[Lambda]},Amat,{bmat,bvec},{cmat,cvec},{x,xa,xb,xmatch}]
CompoundMatrixMethod[{\[Lambda]_/;!NumericQ[\[Lambda]],\[Lambda]\[Lambda]_?NumericQ},Amat_?MatrixQ,bmat_,cmat_,{x_/;!NumericQ[x],xa_?NumericQ,xb_?NumericQ,xmatch_:False},dampfactor_:1]:=CompoundMatrixMethod[{\[Lambda],\[Lambda]\[Lambda]},Amat,{bmat,0},{cmat,0},{x,xa,xb,xmatch}]
CompoundMatrixMethod[{\[Lambda]_/;!NumericQ[\[Lambda]],\[Lambda]\[Lambda]_?NumericQ},Amat_?MatrixQ,{bmat_?MatrixQ,bvec_},{cmat_?MatrixQ,cvec_},{x_/;!NumericQ[x],xa_?NumericQ,xb_?NumericQ,xmatch_:False}]:=Module[{len,subsets,yns,ya,yb,\[Phi]mvec,\[Phi]pvec,bsol,csol,yaaa,ybbb,\[Phi]pa,\[Phi]ma,\[Phi]mb,valsleft,valsright,\[Phi]painit,\[Phi]mbinit,QQ,posint,negint,det,matchpt},
If[(xa<= xmatch<= xb && NumericQ[xmatch]),matchpt=xmatch,matchpt=(xb+xa)/2];
len=Length[Amat];
(*If[!EvenQ[len],Print["Matrix A does not have even dimension"];Abort[]];*)
If[Length[Amat]!=Length[Transpose[Amat]],Print["Matrix A is not a square matrix"];Abort[]];
If[Max@Abs@bvec!=0||Max@Abs@cvec!=0,Print["Warning: BCs are not equal to zero."];];
(*If[!FreeQ[bmat,\[Lambda]],bmat=Extract[#[[2]],Position[#[[1]],a_/;a<0,1]]&@Eigensystem[Amat/.x\[Rule]xa/.\[Lambda]\[Rule]\[Lambda]\[Lambda]]];*)
ya=NullSpace[bmat/.x->xa];
yns=Through[Array[\[FormalY],{len}][x]];
lena=Length[ya];
subsetsa=Subsets[Range[len],{lena}];
(*If[!FreeQ[cmat,\[Lambda]],cmat=Extract[#[[2]],Position[#[[1]],a_/;a>0,1]]&@Eigensystem[Amat/.x\[Rule]xb/.\[Lambda]\[Rule]\[Lambda]\[Lambda]]];*)
yb=NullSpace[cmat/.x->xb];
lenb=Length[yb];
subsetsb=Subsets[Range[len],{lenb}];
If[Length[ya]+Length[yb]!=len,Print["B, C matrices not the right rank"];Abort[]];
\[Phi]mvec=Table[\[FormalPhi]m[i][x],{i,1,Length[subsetsa]}];
\[Phi]pvec=Table[\[FormalPhi]p[i][x],{i,1,Length[subsetsb]}];
bsol=Quiet@Solve[bmat.yns==bvec,yns];
csol=Quiet@Solve[cmat.yns==cvec,yns];
If[Length[bsol]==0||Length[csol]==0,Print["b or c solution failed" ];Return[$Failed]];
yaaa=Transpose[ya+Table[yns,{lena}]/.bsol[[1]]/.Thread[yns->0]];
ybbb=Transpose[yb+Table[yns,{lenb}]/.csol[[1]]/.Thread[yns->0]];
\[Phi]pa=(Det[(yaaa/.x->xa/.\[Lambda]->\[Lambda]\[Lambda])[[#]]]&/@subsetsa);
\[Phi]mb=(Det[(ybbb/.x->xb/.\[Lambda]->\[Lambda]\[Lambda])[[#]]]&/@subsetsb);
valsleft=Select[Eigenvalues[Amat/.x->xa/.\[Lambda]->\[Lambda]\[Lambda]],Re[#]>0&];
valsright=Select[Eigenvalues[Amat/.x->xb/.\[Lambda]->\[Lambda]\[Lambda]],Re[#]<0&];
\[Phi]painit=Thread[Through[Array[\[FormalPhi]p,{Length[subsetsa]}][xa]]==\[Phi]pa];
\[Phi]mbinit=Thread[Through[Array[\[FormalPhi]m,{Length[subsetsb]}][xb]]==\[Phi]mb];
QQa=qMatrix[len,lena]/.\[FormalCapitalA][i_,j_]:>Amat[[i,j]]/.\[Lambda]->\[Lambda]\[Lambda];
QQb=qMatrix[len,lenb]/.\[FormalCapitalA][i_,j_]:>Amat[[i,j]]/.\[Lambda]->\[Lambda]\[Lambda];
posint=NDSolve[{Thread[D[\[Phi]pvec,x]==(QQa-  Total[Re@valsleft] IdentityMatrix[Length[QQa]]).\[Phi]pvec],\[Phi]painit},Array[\[FormalPhi]p,{Length[subsetsa]}],{x,xa,xb}][[1]];
negint=NDSolve[{Thread[D[\[Phi]mvec,x]==(QQb-  Total[Re@valsright] IdentityMatrix[Length[QQb]]).\[Phi]mvec],\[Phi]mbinit},Array[\[FormalPhi]m,{Length[subsetsb]}],{x,xa,xb}][[1]];
det=Total@Table[\[FormalPhi]m[i][x] \[FormalPhi]p[Complement[Range[len],i]][x] (-1)^(Total[Range[lena]+i])//.reprules/.Thread[subsetsa->Range[Length[subsetsa]]]/.Thread[subsetsb->Range[Length[subsetsb]]],{i,subsetsa}];
Exp[-Integrate[Tr[Amat/.\[Lambda]->\[Lambda]\[Lambda]],{x,xa,(xb-xa)/2}]] det/.x->matchpt/.posint/.negint]

ToLinearMatrixForm[eqn_,bcs_?ListQ,depvars_,x_]:=If[bcs=={},ToLinearMatrixForm[eqn,bcs,depvars,{x,0,0}],Print["Please supply the limits for the independent variable"];Return[$Failed]]
ToLinearMatrixForm[eqn_,bcs_?ListQ,depvars_,{x_,xa_,xb_}]:=Module[{yssall,ylen,yssalli,yssi,nODEi,eqns,yns,linearisedeqn,ydvec,fvec,yss,depvars1 ,eqn2,ythr,ythrd,Amat,bcsleft,bcsright,bmat,cmat,bvec,cvec,undifvars,sol,undifsol},
depvars1=Flatten[{depvars}];
If[!ListQ[eqn],eqns={eqn},eqns=eqn];
ylen=Length[depvars1];
undifvars={};
Table[nODEi[i]=Max[0,Cases[eqns,Derivative[m_][depvars1[[i]]][x]:>m,Infinity]];
yssalli[i]=Table[D[depvars1[[i]][x],{x,n}],{n,0,nODEi[i]}];
If[Length[yssalli[i]]==1,AppendTo[undifvars,yssalli[i][[1]]]];
yssi[i]=Most@yssalli[i],
{i,ylen}];
yssall=Flatten[Table[yssalli[i],{i,1,ylen}]];
yss=Flatten[Table[yssi[i],{i,1,ylen}]];
eqn2=eqns/.Thread[yssall->(yssall \[FormalE])];
linearisedeqn=Normal@Series[eqn2,{\[FormalE],0,1}];
If[!(((linearisedeqn/.Equal->Subtract)-(eqn2/.Equal->Subtract) //Simplify)== Table[0,ylen]),(Print["Original Equation not linear, linearized equation: " ,linearisedeqn/.\[FormalE]->1])];
yns=Through[Array[\[FormalY],{Length[yss]}][x]];
ythr=Thread[yss->yns];
ythrd=D[ythr,x];
sol=Solve[Flatten@Join[Select[Thread[D[yss/.ythr,x]==D[yss,x]/.ythr],FreeQ[#,Alternatives@@depvars1]&],{linearisedeqn/.\[FormalE]->1}//.ythr/.ythrd],Join[D[yss/.ythr,x],undifvars]];
If[Length[sol]==0,Print["Error, something has gone wrong!"];Return[$Aborted]];
If[Length[sol]>1,Print["Error, solution space higher dimension than 1"];Return[$Aborted]];
ydvec=D[yss/.ythr,x]/.sol[[1]];
Amat=Coefficient[#,yns]&/@ydvec;
fvec=ydvec/.Thread[yns->0];
If[AnyTrue[fvec,(!PossibleZeroQ[#])&],Print["Inhomogeneous ODE, does not work with CMM at the moment! f = ",fvec]];
If[Length[bcs]!=Length[Amat],Print["Not enough BCs given, returning just the matrix"];Return[Amat]];
If[Length[undifvars]>0,undifsol=(DSolve[Take[sol[[1]],-Length[undifvars]]/.Rule->Equal,depvars,x]//Quiet),undifsol={}];
bcsleft=Select[bcs/.(Thread[yss->yns]/.x->xa)/.(sol[[1]]/.x->xa),!FreeQ[#,\[FormalY]]&];
bcsright=Select[bcs/.(Thread[yss->yns]/.x->xb)/.(sol[[1]]/.x->xb),!FreeQ[#,\[FormalY]]&];
bmat=Coefficient[#,yns/.x->xa]&/@(bcsleft/.Equal->Subtract);
cmat=Coefficient[#,yns/.x->xb]&/@(bcsright/.Equal->Subtract);
bvec=-(bcsleft/.Equal->Subtract)/.Thread[(yns/.x->xa)->0];
cvec=-(bcsright/.Equal->Subtract)/.Thread[(yns/.x->xb)->0];
{Amat,{bmat,bvec},{cmat,cvec},{x,xa,xb}}
]



End[]; (* End Private Context *)

EndPackage[];
