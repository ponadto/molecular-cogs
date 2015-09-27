(* ::Package:: *)

BeginPackage["AffinityPropagation3`"];


myMatrixPlot::usage="";
myMatrixPlotWithNumbers::usage="";
reapForOneAngle::usage="";
affinityPropagation2::usage="";
affinityPropagation2andAHalf::usage="";
packMatrix::usage="";
unpackClustering::usage="";
getMats::usage="";
inputMatrixPath::usage="";
inputdAdKsiPath::usage="";
myTotal::usage="";


PREFIX:="ii-ccl";
INTERACTIONS:={"bonds","angle","tors","vdw","ele"};(*,"oop"};*)
(interactionPosition[INTERACTIONS[[#]]]:=#)&/@Range[Length[INTERACTIONS]];
DIHEDRALS=ToString[NumberForm[#,{4,2}]]&/@Table[i,{i,-179,-60,0.5}];
(dihedralPosition[DIHEDRALS[[#]]]:=#)&/@Range[Length[DIHEDRALS]];
N\[LetterSpace]OF\[LetterSpace]ITERATIONS:=50;
MOL2\[LetterSpace]FILE\[LetterSpace]NAME="";
N\[LetterSpace]OF\[LetterSpace]ATOMS=0;
tmp=StringSplit[#]&/@Import[PREFIX<>".conf","CSV"];
Do[
If[tmp[[i,1,1]]=="mol2FileName",MOL2\[LetterSpace]FILE\[LetterSpace]NAME=tmp[[i,1,2]]];
If[tmp[[i,1,1]]=="nOfAtoms",N\[LetterSpace]OF\[LetterSpace]ATOMS=ToExpression[tmp[[i,1,2]]]];
,{i,1,Length[tmp]}];

atomNumber=1;
tmp=StringSplit[#]&/@Import[MOL2\[LetterSpace]FILE\[LetterSpace]NAME,"CSV"];
ATOMS={};
NUMBERS={};
Do[
If[Length[tmp[[i,1]]]>0&&tmp[[i,1,1]]==ToString[atomNumber],
AppendTo[ATOMS,tmp[[i,1,2]]];
AppendTo[NUMBERS,ToString[atomNumber]];
atomNumber=atomNumber+1;
];
,{i,1,Length[tmp]}];

(*ATOMS:={"N1","H2","H3","H4","C5","C6","I7","I8","H9","H10","H11"};
NUMBERS={"1","2","3","4","5","6","7","8","9","10","11"};*)


Begin["`Private`"];

Needs["SpectralClustering`","spectralClusteringPackage.m"];


myMatrixPlot[mat_,which_,textSize_,size_]:=MatrixPlot[mat,FrameTicks->{{Range[Length[ATOMS[[which]]]],ATOMS[[which]]}\[Transpose],{Range[Length[ATOMS[[which]]]],Rotate[#,90 Degree]&/@ATOMS[[which]]}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Directive[Opacity[1],Bold,Black,textSize-6,FontFamily->"Helvetica"],Mesh->True,ImageSize->size]
myMatrixPlot[mat_]:=myMatrixPlot[mat,Range[Length[ATOMS]],16];
myMatrixPlotWithNumbers[mat_,which_,textSize_,size_]:=MatrixPlot[mat,FrameTicks->{{Range[Length[NUMBERS[[which]]]],NUMBERS[[which]]}\[Transpose],{Range[Length[NUMBERS[[which]]]],Rotate[#,0 Degree]&/@NUMBERS[[which]]}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Directive[Opacity[1],Bold,Black,textSize-6],Mesh->True,ImageSize->size]

inputMatrixPath[angle_,interaction_,what_]:="bootstrapOutput/blockBootstrap"<>what<>"_"<>interaction<>"_"<>ToString[angle]<>"_"<>PREFIX<>".dat";
inputdAdKsiPath[angle_]:="bootstrapOutput/dAdKsiAndZksi="<>angle<>"_"<>PREFIX<>".dat";

reapForOneAngle[angle_,what_]:=Module[{tmp,dim,mats,zksis},

mats=Reap[
Do[
tmp=Import[inputMatrixPath[angle,INTERACTIONS[[i]],what],"Data"];
AppendTo[tmp,{}];
tmp=PadLeft[#,Length[ATOMS]]&/@tmp;
tmp=tmp+tmp\[Transpose];
Sow[tmp];
,
{i,1,Length[INTERACTIONS]}
];
][[2,1]];

mats
];


getMats[interaction_,mats_]:=Module[{positions},

If[MemberQ[INTERACTIONS,interaction],
Return[mats[[interactionPosition[interaction]]]];
];

If[interaction=="non-bonded",
positions=interactionPosition/@{"ele","vdw"};
Return[Total[mats[[positions]]]];
];

If[interaction=="conf",
positions=interactionPosition/@{"bonds","angle","tors"};
Return[Total[mats[[positions]]]];
];

If[interaction=="all",
Return[Total[mats[[All]]]];
];
Print["something went wrong in 'getMat' for interaction = ",interaction];
Assert[1==0];
]

affinityPropagation2[simMat_,dampingCoeff_,howToSetTheDiagonal_]:=Module[{dim,noiseMat,s,a,r,rOld,as,aTmp,diag,aOld,tmp,i,k},
dim=Length[simMat];
(*SeedRandom[1234];*)
noiseMat=Table[RandomReal[]*10^-12*(Max[simMat]-Min[simMat]),{i,1,dim},{j,1,dim}];(*TODO: zastanow sie, czy to cos zmienia*)
s=simMat+(noiseMat+noiseMat\[Transpose])/2;

If[howToSetTheDiagonal==0,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[10^-12*Table[Min[Flatten[(Select[#,#>0&]&/@simMat)]],{dim}]]];
If[howToSetTheDiagonal==1,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[Median/@(Select[#,#!=0&]&/@simMat)]];
If[howToSetTheDiagonal==2,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[Table[Min[Flatten[(Select[#,#>0&]&/@simMat)]],{dim}]]];
If[howToSetTheDiagonal==3,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[Mean/@simMat]];

a=Table[0,{i,1,dim},{j,1,dim}];(*availabilities*)
r=Table[0,{i,1,dim},{j,1,dim}];(*responsibilities*)

Do[
rOld=r;

as=a+s;
(*
r=Table[
s\[LeftDoubleBracket]i,k\[RightDoubleBracket]-Max[as\[LeftDoubleBracket]i,DeleteCases[Range[dim],k]\[RightDoubleBracket]]
,{i,1,dim},{k,1,dim}];
*)
r=s-Map[Max,Table[Drop[#,{k}],{k,1,dim}]&/@as,{2}];
r=(1-dampingCoeff)*r+dampingCoeff*rOld;

aOld=a;
(*
a=Table[
Min[0,r\[LeftDoubleBracket]k,k\[RightDoubleBracket]+Total[Max[0,#]&/@r\[LeftDoubleBracket]All,k\[RightDoubleBracket]]-Max[0,r\[LeftDoubleBracket]k,k\[RightDoubleBracket]]-Max[0,r\[LeftDoubleBracket]i,k\[RightDoubleBracket]]]
,{i,1,dim},{k,1,dim}];
Do[a\[LeftDoubleBracket]k,k\[RightDoubleBracket]=Total[Max[0,#]&/@r\[LeftDoubleBracket]All,k\[RightDoubleBracket]]-Max[0,r\[LeftDoubleBracket]k,k\[RightDoubleBracket]];,{k,1,dim}];
*)
aTmp=Map[Max[#,0]&,r,{2}];
aTmp=aTmp-DiagonalMatrix[Diagonal[aTmp]-Diagonal[r]];
a=Table[Total/@Transpose[aTmp],{k,1,dim}]-aTmp;
diag=Diagonal[a];
a=Map[Min[0,#]&,a,{2}];
a=a-DiagonalMatrix[Diagonal[a]-diag];
(*If[a\[Equal]a2,Print["dobrze jest"],Print["zle jest"]];
Print[MatrixForm[a-a2]];*)
a=(1-dampingCoeff)*a+dampingCoeff*aOld;

,{iter,1,N\[LetterSpace]OF\[LetterSpace]ITERATIONS}];

r+a
];

affinityPropagation2andAHalf[simMat_,dampingCoeff_,howToSetTheDiagonal_]:=Module[{dim,noiseMat,s,a,r,rOld,as,aTmp,diag,aOld,tmp,i,k},
dim=Length[simMat];
(*SeedRandom[1234];*)
noiseMat=Table[RandomReal[]*10^-12*(Max[simMat]-Min[simMat]),{i,1,dim},{j,1,dim}];(*TODO: zastanow sie, czy to cos zmienia*)
s=simMat+(noiseMat+noiseMat\[Transpose])/2;

If[howToSetTheDiagonal==0,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[10^-12*Table[Min[Flatten[(Select[#,#>0&]&/@simMat)]],{dim}]]];
If[howToSetTheDiagonal==1,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[Median/@(Select[#,#!=0&]&/@simMat)]];
If[howToSetTheDiagonal==2,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[Table[Min[Flatten[(Select[#,#>0&]&/@simMat)]],{dim}]]];
If[howToSetTheDiagonal==3,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[Mean/@simMat]];

a=Table[0,{i,1,dim},{j,1,dim}];(*availabilities*)
r=Table[0,{i,1,dim},{j,1,dim}];(*responsibilities*)

Do[
rOld=r;

as=a+s;
(*
r=Table[
s\[LeftDoubleBracket]i,k\[RightDoubleBracket]-Max[as\[LeftDoubleBracket]i,DeleteCases[Range[dim],k]\[RightDoubleBracket]]
,{i,1,dim},{k,1,dim}];
*)
r=s-Map[Max,Table[Drop[#,{k}],{k,1,dim}]&/@as,{2}];
r=(1-dampingCoeff)*r+dampingCoeff*rOld;

aOld=a;
(*
a=Table[
Min[0,r\[LeftDoubleBracket]k,k\[RightDoubleBracket]+Total[Max[0,#]&/@r\[LeftDoubleBracket]All,k\[RightDoubleBracket]]-Max[0,r\[LeftDoubleBracket]k,k\[RightDoubleBracket]]-Max[0,r\[LeftDoubleBracket]i,k\[RightDoubleBracket]]]
,{i,1,dim},{k,1,dim}];
Do[a\[LeftDoubleBracket]k,k\[RightDoubleBracket]=Total[Max[0,#]&/@r\[LeftDoubleBracket]All,k\[RightDoubleBracket]]-Max[0,r\[LeftDoubleBracket]k,k\[RightDoubleBracket]];,{k,1,dim}];
*)
aTmp=Map[Max[#,0]&,r,{2}];
aTmp=aTmp-DiagonalMatrix[Diagonal[aTmp]-Diagonal[r]];
a=Table[Total/@Transpose[aTmp],{k,1,dim}]-aTmp;
diag=Diagonal[a];
a=Map[Min[0,#]&,a,{2}];
a=a-DiagonalMatrix[Diagonal[a]-diag];
(*If[a\[Equal]a2,Print["dobrze jest"],Print["zle jest"]];
Print[MatrixForm[a-a2]];*)
a=(1-dampingCoeff)*a+dampingCoeff*aOld;

,{iter,1,N\[LetterSpace]OF\[LetterSpace]ITERATIONS}];


{r,a}
];

packMatrix[mat_]:=Module[{zeroLine,whichNonZero},
zeroLine=Table[0,{Length[mat]}];
#!=zeroLine&/@mat;
whichNonZero=Flatten[Position[#!=zeroLine&/@mat,True]];
{mat[[whichNonZero,whichNonZero]],whichNonZero}
];
affinityPropagation[simMat_,dampingCoeff_,howToSetTheDiagonal_]:=Module[{dim,noiseMat,s,a,r,rOld,as,aOld,tmp,i,k},
(*
s - similarity matrix;
dampingCoeff - with values between 0 and 1;
howToSetTheDiagonal - expected 1 or 2;
*)
dim=Length[simMat];
(*SeedRandom[1234];*)
noiseMat=Table[RandomReal[]*10^-12*(Max[simMat]-Min[simMat]),{i,1,dim},{j,1,dim}];(*TODO: zastanow sie, czy to cos zmienia*)
s=simMat+(noiseMat+noiseMat\[Transpose])/2;

If[howToSetTheDiagonal==1,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[Median[s[[#,DeleteCases[Range[dim],#]]]]&/@Range[dim]]];
If[howToSetTheDiagonal==2,s=s-DiagonalMatrix[Diagonal[s]]+DiagonalMatrix[Table[1,{dim}]]];

a=Table[0,{i,1,dim},{j,1,dim}];(*availabilities*)
r=Table[0,{i,1,dim},{j,1,dim}];(*responsibilities*)
Do[
rOld=r;
as=a+s;
r=Table[
s[[i,k]]-Max[as[[i,DeleteCases[Range[dim],k]]]]
,{i,1,dim},{k,1,dim}];
r=(1-dampingCoeff)*r+dampingCoeff*rOld;

aOld=a;
a=Table[
Min[
0,
r[[k,k]]+Total[Max[0,#]&/@r[[All,k]]]-Max[0,r[[k,k]]]-Max[0,r[[i,k]]]
],
{i,1,dim},{k,1,dim}];
Do[
a[[k,k]]=Total[Max[0,#]&/@r[[All,k]]]-Max[0,r[[k,k]]];
,{k,1,dim}];
a=(1-dampingCoeff)*a+dampingCoeff*aOld;
,{iter,1,N\[LetterSpace]OF\[LetterSpace]ITERATIONS}];

r+a
];

unpackClustering[clustering_,whichNonZero_,dim_]:=Module[{unpackedClustering},

unpackedClustering=Table[0,{dim}];

Do[
unpackedClustering[[whichNonZero[[i]]]]=clustering[[i]]
,{i,1,Length[clustering]}
];

unpackedClustering
];
myTotal[list_]:=If[Length[list]==0,Return[0],Total[list]];






End[];
EndPackage[];
