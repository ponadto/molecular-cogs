(* ::Package:: *)

BeginPackage["SpectralClustering`"];
largeW::usage="

largeW[numbers_,antinumbers_,simMat_]

INPUT:
numbers - a list of integers, corresponding to nodes (set A of nodes);
antinumbers - list of integers, corresponding to the rest of nodes (set B of nodes);
simMat - list of lists (matrix) of real numbers; simMat\[LeftDoubleBracket]i,j\[RightDoubleBracket] is the weight (similarity) between nodes 'i' and 'j';

OUTPUT:
sum over all vertices joining nodes from 'numbers' with the nodes from 'antinumbers'


";

Ncut::usage="

Ncut[clusters_,simMat_]

INPUT:

OUTPUT:

";
vol::usage="";
plotSimMat::usage="";
plotBlockLikeSimMat::usage="";
getMinChiAndClusters::usage="";
getMinChiAndClusters2::usage=""
getMinChi::usage="";
clusteringBasedOnMembershipMatrix::usage="";
findNumberOfClusters2::usage="";
findCoefficient::usage="";


Begin["`Private`"];
COLORLIST={Blue,Red,Green,Yellow,Cyan,Gray,Orange,Pink,Brown,LightGreen,LightBlue,LightYellow,LightOrange,White};
GRAPHICS="off";
myColorFunction[colorList_]:=#->colorList[[#]] & /@Range[Length[colorList]];
(******************************************************************)
largeW[numbers_,antinumbers_,simMat_]:=Total[Flatten[simMat[[numbers,antinumbers]]]];
vol[numbers_,simMat_]:=largeW[numbers,numbers,simMat];
cut[clusters_,simMat_,k_]:=Module[{numbers,antinumbers},numbers=Flatten[Position[clusters,#]]&/@Range[k];
antinumbers=Complement[Range[Length[simMat]],#]&/@numbers;
0.5*Total[largeW[numbers[[#]],antinumbers[[#]],simMat]&/@Range[k]]
];
Ncut[clusters_,simMat_]:=Module[{numbers,antinumbers,k},
k=Max[clusters];
numbers=Flatten[Position[clusters,#]]&/@Range[k];
If[Length[numbers]==1,
2,
antinumbers=Complement[Range[Length[simMat]],#]&/@numbers;

0.5*Total[(largeW[numbers[[#]],antinumbers[[#]],simMat]/vol[numbers[[#]],simMat])&/@Range[k]]]];
(****************************************************************************************)

plotSimMat[clusters_,simMat_]:=plotSimMat[clusters,simMat,COLORLIST];
plotBlockLikeSimMat[clusters_,simMat_]:=plotBlockLikeSimMat[clusters,simMat,COLORLIST];

plotSimMat[clusters_,simMat_,colorList_]:=Module[{tmp},

tmp=Grid[{
{MatrixPlot[simMat,MaxPlotPoints->Infinity,ImageSize->200,FrameTicks->None]},
{MatrixPlot[Table[clusters[[i]],{j,1,2},{i,1,Length[clusters]}],MaxPlotPoints->Infinity,ColorRules->myColorFunction[colorList],ImageSize->200,FrameTicks->None]}
},
Alignment->Automatic,Frame->True,Spacings->{0,0}];

If[GRAPHICS=="on",
Print[tmp];
];

tmp
];

plotBlockLikeSimMat[clusters_,simMat_,colorList_]:=Module[{map,blockLikeMatrix,sortedClusters},

map=Sort[Range[1,Length[simMat]],clusters[[#1]]<clusters[[#2]]&];
blockLikeMatrix=simMat[[map,map]];
sortedClusters=Sort[clusters];

plotSimMat[sortedClusters,blockLikeMatrix,colorList]
];

distFromHiper[spanningVectors_,point_]:=Module[{orthogonalized,distance},
orthogonalized=Orthogonalize[Flatten[{spanningVectors,{point}},1]];
distance=Dot[point,orthogonalized[[-1]]]*Norm[orthogonalized[[-1]]]];

getMinChiAndClusters[simMat_,k_]:=Module[{matD,matT,eigenvectors,eigenvalues,matY,norms,indicesOfSpanningVectors,restOfTheVectors,distances,nextFarthest,matA,chiMat},
matD=DiagonalMatrix[Total[simMat]];
matT=MatrixPower[matD,-1].simMat;

{eigenvalues,eigenvectors}=Eigensystem[matT];

matY=Transpose[eigenvectors[[1;;k]]];

norms=Norm/@matY;
nextFarthest=Position[norms,Max[norms]][[1]][[1]];
indicesOfSpanningVectors={nextFarthest};
restOfTheVectors=Range[Length[matY]];
Do[
restOfTheVectors=DeleteCases[restOfTheVectors,nextFarthest];
distances=distFromHiper[matY[[indicesOfSpanningVectors]],#]&/@matY[[restOfTheVectors]];
nextFarthest=restOfTheVectors[[Position[distances,Max[distances]][[1]][[1]]]];
AppendTo[indicesOfSpanningVectors,nextFarthest];
,
{i,2,k}
];

matA=MatrixPower[matY[[indicesOfSpanningVectors]],-1];
chiMat=matY.matA;


{Min[chiMat],clusteringBasedOnMembershipMatrix[chiMat],chiMat}
];
getMinChi[simMat_,k_]:=getMinChiAndClusters[simMat,k][[1]];
clusteringBasedOnMembershipMatrix[chiMat_]:=Ordering[#,-1][[1]]&/@chiMat;
findNumberOfClusters2[simMat_,maxClustersTried_,isClusteringValidThreshold_]:=Module[{minChiList,tmp,sorted,nOfClusters,tmpSorted},
If[findCoefficient[simMat]>isClusteringValidThreshold,
If[GRAPHICS=="on",
Print["coefficient = "<>ToString[findCoefficient[simMat]]];
];
Return[{1,{}}];
(*second value is the minChiList, which in this case doesn't make sense*)
];
minChiList=Table[getMinChi[simMat,k],{k,2,maxClustersTried}];
tmp={1+Range[Length[minChiList]],minChiList}\[Transpose];
sorted=Sort[tmp,#1[[2]]>#2[[2]]&];

Assert[sorted[[1,1]]==2];

(*for k=2 clusters minChi should be the closest to 0 (the largest)*)

nOfClusters=If[sorted[[2,2]]>0.5*sorted[[3,2]],sorted[[2,1]],2];

If[GRAPHICS=="on",
Print["minChiList=",minChiList];
Print[ListPlot[tmp,PlotMarkers->{Automatic,Large},AxesLabel->"I chose: "<>ToString[nOfClusters]]];
];

{nOfClusters,minChiList}
];
getMinChiAndClusters2[simMat_,k_]:=Module[{matD,matT,eigenvectors,eigenvalues,matY,norms,indicesOfSpanningVectors,restOfTheVectors,distances,nextFarthest,matA,chiMat},
matD=DiagonalMatrix[Total[simMat]];
matT=MatrixPower[matD,-1].simMat;
{eigenvalues,eigenvectors}=Eigensystem[matT];
matY=eigenvectors[[1;;k]];
matY=matY\[Transpose];
norms=Norm/@matY;
nextFarthest=Position[norms,Max[norms]][[1]][[1]];
indicesOfSpanningVectors={nextFarthest};
restOfTheVectors=Range[Length[matY]];
Do[
restOfTheVectors=DeleteCases[restOfTheVectors,nextFarthest];
distances=distFromHiper[matY[[indicesOfSpanningVectors]],#]&/@matY[[restOfTheVectors]];
nextFarthest=restOfTheVectors[[Position[distances,Max[distances]][[1]][[1]]]];
AppendTo[indicesOfSpanningVectors,nextFarthest];
,
{i,2,k}
];

matA=MatrixPower[matY[[indicesOfSpanningVectors]],-1];
chiMat=matY.matA;

Print[MatrixPlot[matY\[Transpose]]];

{Min[chiMat],ClusteringComponents[matY,k,1],chiMat}
];


(*
to check if a clustering can be made, a 2-clustering is performed, and a 
coefficient = (weights between two clusters)/Min[(weights in cluster 1),(weights in cluster 2)] 
is returned:
*)
findCoefficient[simMat_]:=Module[{clusters,chiMat,numbers,antinumbers},
{clusters,chiMat}=getMinChiAndClusters[simMat,2][[2;;3]];
{numbers,antinumbers}=Flatten[Position[clusters,#]]&/@Range[2];

Print["largeW=",largeW[numbers,antinumbers,simMat]];
Print["vol1=",vol[numbers,simMat]];
Print["vol2=",vol[antinumbers,simMat]];
Print["coefficient=",largeW[numbers,antinumbers,simMat]/Min[vol[numbers,simMat],vol[antinumbers,simMat]]];
largeW[numbers,antinumbers,simMat]/Min[vol[numbers,simMat],vol[antinumbers,simMat]]
];


End[];
EndPackage[];
