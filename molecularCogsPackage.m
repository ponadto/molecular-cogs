(* ::Package:: *)

BeginPackage["MolecularCogs`"];


score::usage="";
separateIntoTwoMatrices::usage="";
adaptSimMat::usage="";
score2::usage="";
assign::usage="";
isSuspectGuilty::usage="";
getClusters::usage="";
assignToClustersBasedOnInteractionMatrix::usage="";
stripExemplars::usage="";
stripAssignments::usage="";
assignToClusters2::usage="";
assignToClusters3::usage="";
assignToClusters4::usage="";
compactClustering::usage="";
mutation::usage="";
crossOver::usage="";
geneticClustering::usage="";
getGeneticClustering::usage="";
calculateScores::usage="";
plotListPlot::usage="";
getIntervalAngles::usage="";
getMatsFromInterval::usage="";
meanForInterval::usage="";
extractPMF::usage="";
extractZksi::usage="";
trapezoidalIntegration::usage="";
extractAccumulatedMat::usage="";
extractGearPMF::usage="";
prepareForErrorListPlot::usage="";
clusteringSimilarity::usage="";
plotListPlotAndStripe::usage="";
plotListPlotAndStripe2::usage="";


SetOptions[ListPlot, 
  TicksStyle -> Directive[FontFamily -> "Helvetica"]];
SetOptions[MatrixPlot, 
  TicksStyle -> Directive[FontFamily -> "Helvetica"]];


Begin["`Private`"];
Needs["AffinityPropagation3`", "affinityPropagationPackage3.m"];
Needs["SpectralClustering`","spectralClusteringPackage.m"];
Needs["ErrorBarPlots`"]


score[clustering_, mat_] :=
  Module[
   {tmp, flatMat, which, plusMat, minusMat, gravityScore, 
    clarityScore, flatMutualMat, totalPlus, totalMinus},
   
   	flatMat = Flatten[mat];
   
   	If[Abs[Total[clustering]] == 
     Length[mat] - 2,(*one gear consists of only one atom*)
    		
    Return[{0, 0}];
    	];
   
   	If[Abs[Total[clustering]] == Length[mat],
    		tmp = Total[flatMat];
    		If[clustering[[1]] > 0,
     			Return[0.5*{tmp/Total[Select[flatMat, # > 0 &]], 0}];
     			,
     			Return[0.5*{tmp/Total[Select[flatMat, # < 0 &]], 0}];
     		];
    	];
   	
   	which = Flatten[Position[clustering, #]] & /@ {-1, 1};
   	plusMat = mat[[which[[2]], which[[2]]]];
   	minusMat = mat[[which[[1]], which[[1]]]];
   totalPlus = myTotal[Select[flatMat, # > 0 &]];
   totalMinus = myTotal[Select[flatMat, # < 0 &]];
   	
   	gravityScore = 
    0.5*(myTotal[myTotal[plusMat]]/totalPlus + 
       myTotal[myTotal[minusMat]]/totalMinus);
   	flatMutualMat = Flatten[mat[[which[[1]], which[[2]]]]];
   	clarityScore = (myTotal[Select[flatMutualMat, # > 0 &]]/
       totalPlus + myTotal[Select[flatMutualMat, # < 0 &]]/totalMinus);
   
   	{gravityScore, -clarityScore}
   ];


separateIntoTwoMatrices[mat_] := Module[{abs},
   abs = Map[Abs, mat, {2}];
   {(mat + abs)/2, (mat - abs)/(-2)}
   ];

adaptSimMat[matrix_] := Module[{min},
   min = Min[Select[Flatten[matrix], # > 0 &]];
   Table[If[matrix[[i, j]] < min, 0.00001*min, matrix[[i, j]]], {i, 1,
      Length[matrix]}, {j, 1, Length[matrix]}]
   ];


score2[clustering_,mat_]:=Module[{tmp,flatMat,which,plusMat,minusMat,lambda,gravityScore,clarityScore,flatMutualMat},

flatMat=Flatten[mat];

(*Print["clustering = ",clustering];*)
If[Abs[Total[clustering]]==Length[mat],
tmp=Total[Total[mat]];
If[tmp>0,
Return[{tmp/Total[Select[flatMat,#>0&]],0}];
,
Return[{tmp/Total[Select[flatMat,#<0&]],0}];
];
];

(*Print[plotBlockLikeSimMat[clustering+5,mat]];*)
which=Flatten[Position[clustering,#]]&/@{-1,1};
plusMat=mat[[which[[2]],which[[2]]]];
minusMat=mat[[which[[1]],which[[1]]]];

lambda=Length[plusMat]/Length[mat];
(*
If[Total[Total[plusMat]]<0||Total[Total[minusMat]]>0,
Return[{0,0}];
];
*)
gravityScore=lambda*Total[Total[plusMat]]/Total[Select[flatMat,#>0&]]+(1-lambda)*Total[Total[minusMat]]/Total[Select[flatMat,#<0&]];
flatMutualMat=Flatten[mat[[which[[1]],which[[2]]]]];
clarityScore=(Total[Select[flatMutualMat,#>0&]]/Total[Select[flatMat,#>0&]]+Total[Select[flatMutualMat,#<0&]]/Total[Select[flatMat,#<0&]]);
(*Print["gravityScore = ",gravityScore];
Print["clarityScore = ",clarityScore];
Print["total score = ",gravityScore-clarityScore];*)

{gravityScore,-clarityScore}
];
assign[results_] := Flatten[Ordering[#1, -1] & /@ results];

isSuspectGuilty[suspect_, assignments_] := Module[{counts},
   counts = Sort[Count[#, suspect] & /@ assignments];
   If[counts[[1]] > 1,
    Return[True];
    ];
   False
   ];

getClusters[assignment_] := 
  Table[Flatten[Position[assignment, node]], {node, 1, 
    Length[assignment]}];

assignToClustersBasedOnInteractionMatrix[interactionMat_, 
   outputFromAP_] := 
  Module[{exemplars, assignments, suspects, absMat, clusters, 
    finalAssignment, clustersSizes, i, tmp, candidates, 
    forFinalResult, finalClustering},
   
   absMat = Map[Abs, interactionMat, {2}];
   
   exemplars = 
    Flatten[Position[Positive[Diagonal[#]], True]] & /@ outputFromAP;
   suspects = Intersection[exemplars[[1]], exemplars[[2]]];
   assignments = assign /@ outputFromAP;
   If[Total[Boole /@ isSuspectGuilty[#, assignments] & /@ suspects] >=
      1,
    Print[MatrixPlot /@ outputFromAP];
    Assert[1 == 0];
    ];
   
   clusters = getClusters /@ assignments;
   Print["clusters = ", clusters];
   clustersSizes = Map[Length, clusters, {2}];
   forFinalResult = 
    Table[{assignments[[#, i]], 
        clustersSizes[[#, assignments[[#, i]]]], 
        absMat[[i, assignments[[#, i]]]]} & /@ Range[2], {i, 1, 
      Length[absMat]}];
   
   (*TODO:*)
   finalClustering = Reap[Do[
       tmp = 
        forFinalResult[[i, 1, 2]]*
         forFinalResult[[i, 2, 
          2]];(*to check if a node belongs to two clusters > 1*)
     
         If[tmp > forFinalResult[[i, 1, 2]] && 
         tmp > forFinalResult[[i, 2, 2]],
        If[
          forFinalResult[[i, 1, 3]] != 0.0 && 
           forFinalResult[[i, 2, 3]] != 0.0 && 
           forFinalResult[[i, 1, 3]] > forFinalResult[[i, 2, 3]], 
          Sow[1], Sow[-1]];
        ,
        If[forFinalResult[[i, 1, 2]] == 1, Sow[-1], Sow[1]];
        ];
       
       , {i, 1, Length[absMat]}]][[2, 1]];
   
   
   finalClustering
   ];

stripExemplars[{exemplarsPlus_, exemplarsMinus_}, {clustersSizesPlus_,
     clustersSizesMinus_}, degrees_] :=
  
  Module[{outputPlus, outputMinus},
   outputPlus = exemplarsPlus;
   outputMinus = exemplarsMinus;
   Do[
    If[clustersSizesPlus[[i]] <= 1 && clustersSizesMinus[[i]] > 1,
     outputPlus = DeleteCases[outputPlus, i];
     ];
    If[clustersSizesMinus[[i]] <= 1 && clustersSizesPlus[[i]] > 1,
     outputMinus = DeleteCases[outputMinus, i];
     ];
    If[clustersSizesPlus[[i]] == 1 && clustersSizesMinus[[i]] == 1,
     If[degrees[[i]] > 0,
       outputMinus = DeleteCases[outputMinus, i];
       ,
       outputPlus = DeleteCases[outputPlus, i];
       ];
     ];
    , {i, 1, Length[clustersSizesPlus]}
    ];
   
   {outputPlus, outputMinus}
   ];

stripAssignments[assignment_, exemplars_] := Module[{output},
   output = Table[0, {Length[assignment]}];
   Do[
    If[MemberQ[exemplars, assignment[[i]]],
      output[[i]] = assignment[[i]];
      ];
    , {i, 1, Length[assignment]}];
   
   output
   ];

assignToClusters2[interactionMat_, outputFromAP_] := 
  Module[{exemplars, assignments, suspects, absMat, degrees, clusters,
     finalAssignment, clustersSizes, i, tmp, candidates, 
    forFinalResult, finalClustering},
   
   absMat = Map[Abs, interactionMat, {2}];
   degrees = Total[interactionMat];
   
   exemplars = 
    Flatten[Position[Positive[Diagonal[#]], True]] & /@ outputFromAP;
   (*Print["exemplars = ",exemplars];*)
   
   suspects = Intersection[exemplars[[1]], exemplars[[2]]];
   assignments = assign /@ outputFromAP;
   (*Print["assignments = ",assignments];*)
   
   If[Total[Boole /@ isSuspectGuilty[#, assignments] & /@ suspects] >=
      1,
    Print[MatrixPlot /@ outputFromAP];
    Print["THERE'S! SOMETHING! WRONG!"];
    Assert[1 == 0];
    ];
   
   clusters = getClusters /@ assignments;
   clustersSizes = Map[Length, clusters, {2}];
   exemplars = stripExemplars[exemplars, clustersSizes, degrees];
   
   assignments = 
    stripAssignments[assignments[[#]], exemplars[[#]]] & /@ Range[2];
   finalClustering = Table[0, {i, 1, Length[interactionMat]}];
   Do[
    If[assignments[[1, i]] == 0,
      finalClustering[[i]] = -1;
      ,
      finalClustering[[i]] = 1;
      ];
    , {i, 1, Length[finalClustering]}];
   
   assignments = Total[assignments];
   
   {finalClustering, assignments}
   ];

assignToClusters3[interactionMat_, outputFromAP_] := Module[
   {exemplars, nonexemplars, assignments, suspects, absMat, degrees, 
    clusters, finalAssignment, clustersSizes, i, tmp, candidates, 
    forFinalResult, finalClustering},
   
   absMat = Map[Abs, interactionMat, {2}];
   degrees = Total[interactionMat];
   
   exemplars = 
    Flatten[Position[Positive[Diagonal[#]], True]] & /@ outputFromAP;
   nonexemplars = 
    Complement[Range[Length[interactionMat]], Flatten[exemplars]];
   suspects = Intersection[exemplars[[1]], exemplars[[2]]];
   assignments = assign /@ outputFromAP;
   
   If[Total[Boole /@ isSuspectGuilty[#, assignments] & /@ suspects] >=
      1,
    Print[MatrixPlot /@ outputFromAP];
    Print["THERE'S! SOMETHING! WRONG!"];
    Assert[1 == 0];
    ];
   
   clusters = getClusters /@ assignments;
   clustersSizes = Map[Length, clusters, {2}];
   exemplars = stripExemplars[exemplars, clustersSizes, degrees];
   
   finalClustering = Table[0, {i, 1, Length[interactionMat]}];
   Do[
    If[assignments[[1, i]] != i && assignments[[2, i]] != i,
      If[interactionMat[[assignments[[1, i]], i]] > 
         Abs[interactionMat[[assignments[[2, i]], i]]],
        finalClustering[[i]] = 1;
        ,
        finalClustering[[i]] = -1;
        ];
      ,
      If[clustersSizes[[1, assignments[[1, i]]]] > 
         clustersSizes[[2, assignments[[2, i]]]],
        finalClustering[[i]] = 1;
        ,
        If[
          clustersSizes[[1, assignments[[1, i]]]] == 
            clustersSizes[[2, assignments[[2, i]]]] && 
           degrees[[i]] > 0,
          finalClustering[[i]] = 1;
          ,
          finalClustering[[i]] = -1;
          ];
        ];
      ];
    , {i, 1, Length[interactionMat]}];
   
   {finalClustering, assignments}
   ];

assignToClusters4[interactionMat_, outputFromAP_] := 
  Module[{exemplars, assignments, suspects, absMat, degrees, clusters,
     finalAssignment, clustersSizes, i, tmp, candidates, 
    forFinalResult, finalClustering},
   
   absMat = Map[Abs, interactionMat, {2}];
   degrees = Total[interactionMat];
   
   exemplars = 
    Flatten[Position[Positive[Diagonal[#]], True]] & /@ outputFromAP;
   suspects = Intersection[exemplars[[1]], exemplars[[2]]];
   assignments = assign /@ outputFromAP;
   If[Total[Boole /@ isSuspectGuilty[#, assignments] & /@ suspects] >=
      1,
    Print[MatrixPlot /@ outputFromAP];
    Print["THERE'S! SOMETHING! WRONG!"];
    Assert[1 == 0];
    ];
   
   finalClustering = Table[0, {i, 1, Length[interactionMat]}];
   Do[
    If[outputFromAP[[1, i, assignments[[1, i]]]] > 
       outputFromAP[[2, i, assignments[[2, i]]]],
      finalClustering[[i]] = 1;
      ,
      finalClustering[[i]] = -1;
      ];
    , {i, 1, Length[interactionMat]}];
   
   
   {finalClustering, assignments}
   ];

compactClustering[mat_, howToSetTheDiagonal_, clusterAssignment_, 
   coefficient_] :=
  Module[
   {packedMat, which, simMats, outputFromAP2, finalClustering, scores},
   
   	{packedMat, which} = packMatrix[mat];
   
   	simMats = {packedMat, -packedMat};
   
   	outputFromAP2 = 
    affinityPropagation2[#, coefficient, howToSetTheDiagonal] & /@ 
     simMats;
   	finalClustering = clusterAssignment[packedMat, outputFromAP2][[1]];
   
   	{finalClustering, packedMat, which}
   ];
mutation[clustering_] := Module[{tmp, position},
  tmp = clustering;
  position = RandomInteger[{1, Length[clustering]}];
  tmp[[position]] = -clustering[[position]];
  tmp
  ]
crossOver[{pater_, mater_}, population_] := 
  Module[{paterCopy, materCopy, positions},
   paterCopy = population[[pater]];
   materCopy = population[[mater]];
   
   positions = 
    RandomSample[Range[Length[paterCopy]], Floor[Length[paterCopy]/2]];
   
   materCopy[[positions]] = population[[pater]][[positions]];
   
   materCopy
   ];
geneticClustering[packedMat_, convergenceRequirement_, 
   populationSize_] := 
  Module[{bestClustering, clusterAssignments, bestScore, bestScores, 
    simMats, population, parents, offspring, scores, outputFromAP2, 
    nOfIterationsWithNoImprovement, whichMax, mutants},
   
   simMats = {packedMat, -packedMat};
   clusterAssignments = {assignToClusters3, assignToClusters4};
   bestScore = -1;
   bestClustering = 0;
   nOfIterationsWithNoImprovement = 0;
   population = Reap[
      Do[
       outputFromAP2 = 
        affinityPropagation2[#, 0.1*coefficient, 
           howToSetTheDiagonal] & /@ simMats;
       Sow[clusterAssignments[[i]][packedMat, outputFromAP2][[1]]];
       ,
       {i, 1, Length[clusterAssignments]}, {howToSetTheDiagonal, 0, 
        3}, {coefficient, 1, 5}]
      ][[2, 1]];
   population = 
    Union[population, Table[#, {Length[packedMat]}] & /@ {-1, 1}];
   
   bestScores = {};
   While[nOfIterationsWithNoImprovement < convergenceRequirement,
    scores = Total[score[#, packedMat]] & /@ population;
    parents = 
     RandomSample[(Exp[2*#] & /@ scores) -> Range[Length[population]],
         2] & /@ Range[Length[population]];
    offspring = crossOver[#, population] & /@ parents;
    population = Union[population, offspring];
    mutants = 
     mutation /@ 
      RandomChoice[population, Floor[0.01*Length[population]]];
    population = Union[population, mutants];
    
    scores = Total[score[#, packedMat]] & /@ population;
    whichMax = Ordering[scores, -1][[1]];
    If[scores[[whichMax]] > bestScore,
     bestScore = scores[[whichMax]];
     bestClustering = population[[whichMax]];
     nOfIterationsWithNoImprovement = 0;
     ,
     nOfIterationsWithNoImprovement += 1;
     ];
    AppendTo[bestScores, bestScore];
    
    If[Length[population] > populationSize,
     population = 
       RandomSample[(N[Exp[#]]&/@(Max[{-2,#}]&/@scores)) -> population, 
        populationSize];
     ,
     population = 
       Join[population, 
        Table[RandomChoice[{-1, 1}, 
          Length[population[[1]]]], {populationSize - 
           Length[population]}]];
     ];
    ];
   
   {bestScore, bestClustering}
   ];
getGeneticClustering[angle_, interaction_, 
   nOfGeneticIterationsWithNoImprovement_, populationSize_] :=
  
  Module[
   {mats, mat, packedMat, which},
   mats = reapForOneAngle[angle, "Mean"];
   mat = mats[[interactionPosition[interaction]]];
   	{packedMat, which} = packMatrix[mat];
   geneticClustering[packedMat, nOfGeneticIterationsWithNoImprovement,
     populationSize]
   ];
calculateScores[dihedral_, interaction_] :=
  Module[
   {mats, particularInteractionMat, finalClustering, myResults, 
    scores, packedMat, which, fullMat, whichLeft, whichRight, 
    pullLeftFromGear, pullRightFromGear, interGearPullLeft, 
    interGearPullRight, pullLeftParticularTotal, 
    pullRightParticularTotal, pullLeftFullTotal, pullRightFullTotal},
  PrintTemporary[dihedral];
   {particularInteractionMat, fullMat} = 
    getMats[#, reapForOneAngle[dihedral, "Mean"]] & /@ {interaction, 
      "all"};
   
   {packedMat, which} = packMatrix[particularInteractionMat];
   
   myResults = geneticClustering[packedMat, 10, 200];
   
   scores = {score[#, packedMat], #} & /@ 
     Tuples[{-1, 1}, Length[packedMat]];
   scores = Sort[scores, Total[#1[[1]]] > Total[#2[[1]]] &];
   
   finalClustering = 
    unpackClustering[myResults[[2]], which, 
     Length[particularInteractionMat]];
   
   {whichLeft, whichRight} = 
    Flatten[Position[finalClustering, #]] & /@ {-1, 1};
   {pullLeftFromGear, pullRightFromGear} = 
    0.5*myTotal[
        myTotal[particularInteractionMat[[#, #]]]] & /@ {whichLeft, 
      whichRight};
   {interGearPullLeft, interGearPullRight} = 
    myTotal /@ {Select[
       Flatten[particularInteractionMat[[whichLeft, whichRight]]], # <
          0 &], Select[
       Flatten[particularInteractionMat[[whichLeft, whichRight]]], # <
          0 &]};
   {pullLeftParticularTotal, pullRightParticularTotal} = 
    0.5*{myTotal[Select[Flatten[particularInteractionMat], # < 0 &]], 
      myTotal[Select[Flatten[particularInteractionMat], # > 0 &]]};
   {pullLeftFullTotal, pullRightFullTotal} = 
    0.5*{Total[Select[Flatten[fullMat], # < 0 &]], 
      Total[Select[Flatten[fullMat], # > 0 &]]};
   
   {{myResults[[1]], Total[scores[[1, 1]]]}, 
    finalClustering, {pullLeftFromGear/pullLeftParticularTotal, 
     pullLeftFromGear/pullLeftFullTotal, 
     interGearPullLeft/pullLeftParticularTotal, 
     interGearPullLeft/pullLeftFullTotal}, {pullRightFromGear/
      pullRightParticularTotal, pullRightFromGear/pullRightFullTotal, 
     interGearPullRight/pullRightParticularTotal, 
     interGearPullRight/pullRightFullTotal}}
   ];
plotListPlot[textSize_, wyniki_, interaction_, size_] := Module[{},
  ListPlot[{{ToExpression /@ DIHEDRALS, 
      wyniki[[All, 1]]}\[Transpose], {ToExpression /@ DIHEDRALS, 
      wyniki[[All, 2]]}\[Transpose]}, Filling -> Axis, 
   PlotRange -> {All, {-0.8, 1}}, Joined -> True, 
   PlotLegends -> 
    Placed[PointLegend[Automatic, 
      Style[#, Bold, textSize - 6, 
         FontFamily -> "Helvetica"] & /@ {"genetic clustering", 
        "best clustering"}(*,LegendFunction\[Rule](Framed[#,
      RoundingRadius\[Rule]1,FrameStyle\[Rule]LightGray]&)*)], {0.4, 
      0.2}], PlotMarkers -> {Automatic, Tiny}, Frame -> True, 
   FrameLabel -> {Style["\[Xi] [degrees]", Bold, textSize, 
      FontFamily -> "Helvetica"], 
     Style["SCORE", textSize, Bold, FontFamily -> "Helvetica"]}, 
   LabelStyle -> 
    Directive[Black, textSize - 6, Bold, FontFamily -> "Helvetica"], 
   ImageSize -> size]
  ]
getIntervalAngles[max_, min_] := Module[{i},
  Table[ToString[NumberForm[i, {100, 2}]], {i, ToExpression[min], 
     ToExpression[max], 0.5}][[;; -1]]
  ]
getMatsFromInterval[max_, min_, interaction_] := 
 Module[{angles, mats, varMats, diff},
  angles = getIntervalAngles[max, min];
  mats = (getMats[interaction, reapForOneAngle[#, "Mean"]]) & /@ 
    angles;
  varMats = (getMats[interaction, reapForOneAngle[#, "Variance"]]) & /@
     angles;
  diff = Abs[
    ToExpression[angles[[2]]] - 
     ToExpression[
      angles[[1]]]];(*I'm assuming that all diffs are the same*)
  
  {angles, mats, varMats, diff}
  ]
meanForInterval[mats_, varMats_, diff_] := 
 Module[{mat, matSum, varMat},
  mat = 0.5*diff*(mats[[1]] + 2*Total[mats[[2 ;; -2]]] + mats[[-1]]);
  varMat = 
   diff*diff*(0.25*varMats[[1]] + Total[varMats[[2 ;; -2]]] + 
      0.25*varMats[[-1]]);
  {mat, varMat}
  ]
extractPMF[max_, min_, which_] := Module[{angles, lambdaFunction},
   angles = getIntervalAngles[max, min];
   lambdaFunction[angle_] := Module[{tmp},
     tmp = 
      ToExpression /@ (Import[inputdAdKsiPath[angle]][[
         which, {3, 5}]]);
     {tmp[[1]], tmp[[2]]}
     ];
   lambdaFunction /@ angles
   ];
extractZksi[max_, min_] := Module[{angles, lambdaFunction},
   angles = getIntervalAngles[max, min];
   lambdaFunction[angle_] := Module[{tmp},
     tmp = 
      ToExpression /@ (Import[inputdAdKsiPath[angle]][[2, {3, 5}]]);
     {tmp[[1]], tmp[[2]]}
     ];
   lambdaFunction /@ angles
   ];
trapezoidalIntegration[{array_, var_}, diff_] := Module[{a, errorA},
   a = {0};
   errorA = {0};
   Do[
    AppendTo[a, a[[-1]] + 0.5*(array[[i]] + array[[i + 1]])*diff];
    AppendTo[errorA, 
     errorA[[-1]] + 0.25*diff*diff*(4*var[[i]] + 4*var[[i + 1]])];
    , {i, 1, Length[array] - 1}];
   errorA[[-1]] = errorA[[-1]] - 0.25*diff*diff*3*var[[-1]];
   {a, errorA}
   ];
extractAccumulatedMat[max_, min_, interaction_] := 
  Module[{angles, mats, varMats, diff, accumulatedMat, rewoundedMats, 
    i, j, k},
   {angles, mats, varMats, diff} = 
    getMatsFromInterval[max, min, interaction];
   accumulatedMat = 
    Table[trapezoidalIntegration[{mats[[All, i, j]], 
       varMats[[All, i, j]]}, 0.5], {i, 1, Length[mats[[1]]]}, {j, 1, 
      Length[mats[[1]]]}];
   rewoundedMats = 
    Table[Table[#[[k, i, j]], {k, 1, Length[angles]}] & /@ {mats, 
       varMats}, {i, 1, Length[mats[[1]]]}, {j, 1, Length[mats[[1]]]}];
   {rewoundedMats, accumulatedMat}
   ];
extractGearPMF[accumulatedMat_, forwardGear_, reverseGear_] := 
  Module[{forwardPull, reversePull, gearGrinding, totalPull, i, 
    dimension},
   dimension = Dimensions[accumulatedMat][[4]];
   forwardPull = 
    0.5*Table[
        myTotal[myTotal[
          accumulatedMat[[forwardGear, forwardGear, #, i]]]], {i, 1, 
         dimension}] & /@ {1, 2};
   reversePull = 
    0.5*Table[
        myTotal[myTotal[
          accumulatedMat[[reverseGear, reverseGear, #, i]]]], {i, 1, 
         dimension}] & /@ {1, 2};
   gearGrinding = 
    Table[myTotal[
        myTotal[Map[Abs, 
          accumulatedMat[[forwardGear, reverseGear, #, 
           i]], {2}]]], {i, 1, dimension}] & /@ {1, 2};
   totalPull = 
    0.5*Table[
        myTotal[myTotal[accumulatedMat[[All, All, #, i]]]], {i, 1, 
         dimension}] & /@ {1, 2};
   {forwardPull, reversePull, gearGrinding, totalPull}
   ];
prepareForErrorListPlot[integral_] := Module[{angles, i},
  angles = Range[-179, -60, 0.5];
  Table[{{angles[[i]], integral[[1, i]]}, 
    ErrorBar[Sqrt[integral[[2, i]]]]}, {i, 1, Length[angles]}]
  ]
clusteringSimilarity[p_, q_, nOfAtoms_] := N[Dot[p, q]/nOfAtoms];
(*SetOptions[ListPlot, 
  TicksStyle -> Directive[FontFamily -> "Helvetica"]];
SetOptions[MatrixPlot, 
  TicksStyle -> Directive[FontFamily -> "Helvetica"]];*)
plotListPlotAndStripe[interaction_, prefix_, textSize_, size_] := 
 Module[{fileName, wyniki, angles, mats, varMats, diff, max, min, p1, 
   p2, mat, varMat, packedMat, which, myResults, p3, accumulatedMat, 
   unpackedClustering, reverseGear, forwardGear, p4, p5, forwardPull, 
   reversePull, gearGrinding, totalPull, upperLimit, lowerLimit, step,
    ticks, p6, p7, p8, map, zeroLeftPosition, zeroRightPosition},
  fileName = "wyniki_" <> interaction <> "_"<>PREFIX<>".str";
  If[FileExistsQ[fileName],
Print["Reading ",fileName];
   wyniki = 
     Uncompress[
      Import[fileName, "String"]];
   ,
   wyniki = (calculateScores[#, interaction]) & /@ DIHEDRALS;
   \!\(TraditionalForm\`\(Export[fileName, Compress[wyniki], "\<String\>"];\)\)
   ];

Print[interaction];
  
 p1=plotListPlot[textSize+3.5,wyniki[[All,1]],interaction,size*2.9/2.09];
Export[prefix<>"Scores_"<>PREFIX<>".pdf",p1];
Print[p1];
(*ListPlot[{wyniki\[LeftDoubleBracket]All,3,1\[RightDoubleBracket],wyniki\[LeftDoubleBracket]All,4,1\[RightDoubleBracket],wyniki\[LeftDoubleBracket]All,3,3\[RightDoubleBracket],wyniki\[LeftDoubleBracket]All,4,3\[RightDoubleBracket]},Joined\[Rule]True,PlotStyle\[Rule]{Directive[Blue,Thick],Directive[Orange,Thick],Directive[Blue,Thick,Dashed],Directive[Orange,Thick,Dashed]},PlotRange\[Rule]{All,{0,1.1}}]
ListPlot[{wyniki\[LeftDoubleBracket]All,3,2\[RightDoubleBracket],wyniki\[LeftDoubleBracket]All,4,2\[RightDoubleBracket],wyniki\[LeftDoubleBracket]All,3,4\[RightDoubleBracket],wyniki\[LeftDoubleBracket]All,4,4\[RightDoubleBracket]},Joined\[Rule]True,PlotStyle\[Rule]{Directive[Blue,Thick],Directive[Orange,Thick],Directive[Blue,Thick,Dashed],Directive[Orange,Thick,Dashed]},PlotRange\[Rule]{All,{0,1.1}}]*)
p2=Labeled[MatrixPlot[wyniki[[All,2]]\[Transpose],MaxPlotPoints->Infinity,FrameTicks->{{Range[Length[wyniki[[All,2]]\[Transpose]]],Rotate[#,0 Degree]&/@ATOMS}\[Transpose],{{1,234,120,160,200,80,40},{-180,-60,-120,-100,-80,-140,-160}}\[Transpose]},PlotRangePadding->0,LabelStyle->Directive[Black,textSize-6,Bold,FontFamily->"Helvetica"],AspectRatio->0.5,Mesh->{True,False}],Style["\[Xi] [degrees]",Bold,textSize,FontFamily->"Helvetica"],ImageSize->size*2.9/2.09];
Export[prefix<>"Stripes_"<>PREFIX<>".pdf",p2];
Print[p2];

max="-62.50";
min="-172.50";
{angles,mats,varMats,diff}=getMatsFromInterval[max,min,interaction];
{mat,varMat}=meanForInterval[mats,varMats,diff];


{packedMat,which}=packMatrix[mat];
myResults=geneticClustering[packedMat,10,200];

(*{whichLeft,whichRight}=Flatten[Position[myResults\[LeftDoubleBracket]2\[RightDoubleBracket],#]]&/@{1,-1}
0.5*Total[Total[packedMat\[LeftDoubleBracket]whichLeft,whichLeft\[RightDoubleBracket]]]
0.5*Total[Total[packedMat\[LeftDoubleBracket]whichRight,whichRight\[RightDoubleBracket]]]
Total[Total[packedMat\[LeftDoubleBracket]whichLeft,whichRight\[RightDoubleBracket]]]
plotBlockLikeSimMat[myResults\[LeftDoubleBracket]2\[RightDoubleBracket],packedMat]
score[myResults\[LeftDoubleBracket]2\[RightDoubleBracket],myResults\[LeftDoubleBracket]2\[RightDoubleBracket]]*)
p3=myMatrixPlot[mat,Range[Length[mat]],textSize,size];
Export[prefix<>"IntMat_"<>PREFIX<>".pdf",p3];
Print[p3];
(*
p4=MatrixPlot[mat\[LeftDoubleBracket]elements,elements\[RightDoubleBracket],FrameTicks->{{Range[Length[ATOMS[[which]]]],ATOMS[[which]]}\[Transpose],{Range[Length[ATOMS[[which]]]],ATOMS[[which]]}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Opacity[1],Mesh\[Rule]All,GridLines\[Rule]{{3,7},{4,8}},GridLinesStyle\[Rule]{Directive[Thick,Black],Directive[Thick,Black]},Method\[Rule]{"GridLinesInFront"\[Rule]True}];
*)

{mats,accumulatedMat}=extractAccumulatedMat[max,min,interaction];
mat=accumulatedMat[[All,All,1,-1]];
{packedMat,which}=packMatrix[mat];
myResults=geneticClustering[packedMat,10,200];
unpackedClustering=unpackClustering[myResults[[2]],which,Length[mat]];
{reverseGear,forwardGear}=Flatten[Position[unpackedClustering,#]]&/@{-1,1};
p5=Labeled[MatrixPlot[{unpackedClustering},FrameTicks->{None,{Range[Length[ATOMS]],ATOMS}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Directive[Opacity[1],Bold,Black,textSize-6,FontFamily->"Helvetica"],Mesh->True],Style["SCORE = "<>ToString[NumberForm[score[myResults[[2]],packedMat][[1]],{4,3}]],FontFamily->"Helvetica",FontSize->textSize,Bold],Top,ImageSize->size*2.9/2.09];
Export[prefix<>"Stripe_"<>PREFIX<>".pdf",p5];
Print[p5];

map=Sort[Range[1,Length[mat]],unpackedClustering[[#1]]<unpackedClustering[[#2]]&];
p4=If[Length[Flatten[Position[unpackedClustering,0]]]==0,
zeroLeftPosition=Min[Flatten[Position[unpackedClustering[[map]],1]]];
MatrixPlot[mat[[map,map]],FrameTicks->{{Range[Length[ATOMS[[map]]]],ATOMS[[map]]}\[Transpose],{Range[Length[ATOMS[[map]]]],Rotate[#,90 Degree]&/@ATOMS[[map]]}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Directive[Opacity[1],Bold,Black,textSize-6,FontFamily->"Helvetica"],Mesh->All,GridLines->{{zeroLeftPosition-1},{Length[ATOMS]-zeroLeftPosition+1}},GridLinesStyle->{Directive[Thick,Black],Directive[Thick,Black]},Method->{"GridLinesInFront"->True},ImageSize->size]
,
{zeroLeftPosition,zeroRightPosition}=#[Flatten[Position[unpackedClustering,0]]]&/@{Min,Max};
If[prefix=="ele"||prefix=="vdw",
{zeroLeftPosition,zeroRightPosition}={zeroLeftPosition,zeroRightPosition}+1;
];
MatrixPlot[mat[[map,map]],FrameTicks->{{Range[Length[ATOMS[[map]]]],ATOMS[[map]]}\[Transpose],{Range[Length[ATOMS[[map]]]],Rotate[#,90 Degree]&/@ATOMS[[map]]}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Directive[Opacity[1],Bold,Black,textSize-6,FontFamily->"Helvetica"],Mesh->All,GridLines->{{zeroLeftPosition-1,zeroRightPosition},{Length[ATOMS]-zeroLeftPosition+1,Length[ATOMS]-zeroRightPosition}},GridLinesStyle->{Directive[Thick,Black],Directive[Thick,Black]},Method->{"GridLinesInFront"->True},ImageSize->size]
];
If[prefix=="total"||prefix=="conf",
p4=p3;
];
Export[prefix<>"IntMatClustered_"<>PREFIX<>".pdf",p4];
Print[p4];

p6=ListPlot[{ToExpression/@DIHEDRALS,Table[clusteringSimilarity[unpackedClustering,wyniki[[All,2]][[i]],Length[ATOMS]],{i,1,Length[wyniki[[All,2]]]}]}\[Transpose],Joined->True,Frame->True,FrameLabel->{Style["\[Xi] [degrees]",Bold,textSize,FontFamily->"Helvetica"],Style["SIMILARITY",textSize,Bold,FontFamily->"Helvetica"]},LabelStyle->Directive[Black,textSize-6,Bold,FontFamily->"Helvetica"],PlotStyle->Thick,PlotRange->{All,{-1,1}},Filling->Bottom,ImageSize->size];
Export[prefix<>"Similarity_"<>PREFIX<>".pdf",p6];
Print[p6];

max="-60.00";
min="-179.00";
{mats,accumulatedMat}=extractAccumulatedMat[max,min,interaction];
(*{forwardPull,reversePull,gearGrinding,totalPull}=extractGearPMF[mats,forwardGear,reverseGear];
upperLimit=Max[Flatten[{forwardPull\[LeftDoubleBracket]1,All\[RightDoubleBracket],reversePull\[LeftDoubleBracket]1,All\[RightDoubleBracket],gearGrinding\[LeftDoubleBracket]1,All\[RightDoubleBracket]}]];
lowerLimit=Min[Flatten[{forwardPull\[LeftDoubleBracket]1,All\[RightDoubleBracket],reversePull\[LeftDoubleBracket]1,All\[RightDoubleBracket],gearGrinding\[LeftDoubleBracket]1,All\[RightDoubleBracket]}]];
step=(upperLimit-lowerLimit)/4;
ticks=Range[lowerLimit-step,upperLimit+step,step];
p7=ListPlot[{{Range[-179,-60,0.5],forwardPull\[LeftDoubleBracket]1,All\[RightDoubleBracket]}\[Transpose],{Range[-179,-60,0.5],reversePull\[LeftDoubleBracket]1,All\[RightDoubleBracket]}\[Transpose],{Range[-179,-60,0.5],gearGrinding\[LeftDoubleBracket]1,All\[RightDoubleBracket]}\[Transpose],{Range[-179,-60,0.5],totalPull\[LeftDoubleBracket]1,All\[RightDoubleBracket]}\[Transpose]},PlotLegends\[Rule]Placed[PointLegend[{Style["FGC",Bold,textSize-6],Style["RGC",Bold,textSize-6],Style["GG",Bold,textSize-6],Style["total",Bold,textSize-6]},LegendMarkerSize\[Rule]textSize-6,LegendMarkers\[Rule]"Filled"],Above],Joined\[Rule]True,PlotRange\[Rule]All,PlotStyle\[Rule]Thick,Frame\[Rule]{{True,False},{True,False}},FrameLabel\[Rule]{Style["\[Xi] [degrees]",Bold,textSize],Style["kcal/mol",textSize,Bold]},FrameTicks\[Rule]{{Automatic,None},{{-180,-160,-140,-120,-100,-80,-60},None}},LabelStyle\[Rule]Directive[Black,textSize-6,Bold],AspectRatio\[Rule]0.8,ImageSize\[Rule]size,TicksStyle\[Rule]Directive[textSize-6,Bold,FontFamily\[Rule]"Helvetica"]];
Export[prefix<>"PMFderivative.pdf",p7,"PDF"];
Print[p7];*)

{forwardPull,reversePull,gearGrinding,totalPull}=extractGearPMF[accumulatedMat,forwardGear,reverseGear];
p8=ErrorListPlot[prepareForErrorListPlot/@{totalPull,gearGrinding,forwardPull,reversePull},PlotLegends->Placed[PointLegend[{Style["total",Bold,textSize-6,FontFamily->"Helvetica"],Style["gear\ngrinding",Bold,textSize-6,FontFamily->"Helvetica"],Style["reverse cogs\ncontribution",Bold,textSize-6,FontFamily->"Helvetica"],Style["forward cogs\ncontribution",Bold,textSize-6,FontFamily->"Helvetica"]},LegendMarkerSize->textSize-6,LegendMarkers->"Filled"],Right],Joined->True,PlotRange->All,Frame->{{True,False},{True,False}},FrameLabel->{Style["\[Xi] [degrees]",Bold,textSize,FontFamily->"Helvetica"],Style["kcal/mol",textSize,Bold]},PlotStyle->{Darker[Green],Purple,Orange,Blue},FrameTicks->{{Automatic,None},{{-180,-160,-140,-120,-100,-80,-60},None}},LabelStyle->Directive[Black,textSize-6,Bold,FontFamily->"Helvetica"],AspectRatio->0.8,ImageSize->size,TicksStyle->Directive[textSize-6,Bold,FontFamily->"Helvetica"]];
Export[prefix<>"PMF_"<>PREFIX<>".pdf",p8];
Print[p8];

 fileName = "forFigures_" <> interaction <> "_"<>PREFIX<>".str";
\!\(TraditionalForm\`\(Export[fileName, Compress[{wyniki, mat, unpackedClustering, myResults, packedMat, map, totalPull, gearGrinding, forwardPull, reversePull}], "\<String\>"];\)\)

{p1,p2,p3,p4,p5,p6,p8}
  ];











plotListPlotAndStripe2[interaction_, prefix_, textSize_, size_] := 
 Module[{fileName, wyniki, angles, mats, varMats, diff, max, min, p1, 
   p2, mat, varMat, packedMat, which, myResults, p3, accumulatedMat, 
   unpackedClustering, reverseGear, forwardGear, p4, p5, forwardPull, 
   reversePull, gearGrinding, totalPull, upperLimit, lowerLimit, step,
    ticks, p6, p7, p8, map, zeroLeftPosition, zeroRightPosition},

fileName = "forFigures_" <> interaction <> "_"<>PREFIX<>".str"; 
Print["Reading ",fileName];
{wyniki,mat,unpackedClustering,myResults,packedMat,map,totalPull,gearGrinding,forwardPull,reversePull}=Uncompress[Import[fileName, "String"]];
 
  
 p1=plotListPlot[textSize+3.5,wyniki[[All,1]],interaction,size*2.9/2.09];
Export[prefix<>"Scores.pdf",p1];
Print[p1];

p2=Labeled[MatrixPlot[wyniki[[All,2]]\[Transpose],MaxPlotPoints->Infinity,FrameTicks->{{Range[Length[wyniki[[All,2]]\[Transpose]]],Rotate[#,0 Degree]&/@ATOMS}\[Transpose],{{1,234,120,160,200,80,40},{-180,-60,-120,-100,-80,-140,-160}}\[Transpose]},PlotRangePadding->0,LabelStyle->Directive[Black,textSize-6,Bold,FontFamily->"Helvetica"],AspectRatio->0.5,Mesh->{True,False}],Style["\[Xi] [degrees]",Bold,textSize,FontFamily->"Helvetica"],ImageSize->size*2.9/2.09];
Export[prefix<>"Stripes.pdf",p2];
Print[p2];

p3=myMatrixPlot[mat,Range[Length[mat]],textSize,size];
Export[prefix<>"IntMat.pdf",p3];
Print[p3];

p5=Labeled[MatrixPlot[{unpackedClustering},FrameTicks->{None,{Range[Length[ATOMS]],ATOMS}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Directive[Opacity[1],Bold,Black,textSize-6,FontFamily->"Helvetica"],Mesh->True],Style["SCORE = "<>ToString[NumberForm[score[myResults[[2]],packedMat][[1]],{4,3}]],FontFamily->"Helvetica",FontSize->textSize,Bold],Top,ImageSize->size*2.9/2.09];
Export[prefix<>"Stripe.pdf",p5];
Print[p5];

map=Sort[Range[1,Length[mat]],unpackedClustering[[#1]]<unpackedClustering[[#2]]&];
p4=If[Length[Flatten[Position[unpackedClustering,0]]]==0,
zeroLeftPosition=Min[Flatten[Position[unpackedClustering[[map]],1]]];
MatrixPlot[mat[[map,map]],FrameTicks->{{Range[Length[ATOMS[[map]]]],ATOMS[[map]]}\[Transpose],{Range[Length[ATOMS[[map]]]],Rotate[#,90 Degree]&/@ATOMS[[map]]}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Directive[Opacity[1],Bold,Black,textSize-6,FontFamily->"Helvetica"],Mesh->All,GridLines->{{zeroLeftPosition-1},{Length[ATOMS]-zeroLeftPosition+1}},GridLinesStyle->{Directive[Thick,Black],Directive[Thick,Black]},Method->{"GridLinesInFront"->True},ImageSize->size]
,
{zeroLeftPosition,zeroRightPosition}=#[Flatten[Position[unpackedClustering,0]]]&/@{Min,Max};
If[prefix=="ele"||prefix=="vdw",
{zeroLeftPosition,zeroRightPosition}={zeroLeftPosition,zeroRightPosition}+1;
];
MatrixPlot[mat[[map,map]],FrameTicks->{{Range[Length[ATOMS[[map]]]],ATOMS[[map]]}\[Transpose],{Range[Length[ATOMS[[map]]]],Rotate[#,90 Degree]&/@ATOMS[[map]]}\[Transpose]},FrameStyle->Opacity[0],FrameTicksStyle->Directive[Opacity[1],Bold,Black,textSize-6,FontFamily->"Helvetica"],Mesh->All,GridLines->{{zeroLeftPosition-1,zeroRightPosition},{Length[ATOMS]-zeroLeftPosition+1,Length[ATOMS]-zeroRightPosition}},GridLinesStyle->{Directive[Thick,Black],Directive[Thick,Black]},Method->{"GridLinesInFront"->True},ImageSize->size]
];
If[prefix=="total"||prefix=="conf",
p4=p3;
];
Export[prefix<>"IntMatClustered.pdf",p4];
Print[p4];

p6=ListPlot[{ToExpression/@DIHEDRALS,Table[clusteringSimilarity[unpackedClustering,wyniki[[All,2]][[i]],Length[ATOMS]],{i,1,Length[wyniki[[All,2]]]}]}\[Transpose],Joined->True,Frame->True,FrameLabel->{Style["\[Xi] [degrees]",Bold,textSize,FontFamily->"Helvetica"],Style["SIMILARITY",textSize,Bold,FontFamily->"Helvetica"]},LabelStyle->Directive[Black,textSize-6,Bold,FontFamily->"Helvetica"],PlotStyle->Thick,PlotRange->{All,{-1,1}},Filling->Bottom,ImageSize->size];
Export[prefix<>"Similarity.pdf",p6];
Print[p6];

p8=ErrorListPlot[prepareForErrorListPlot/@{totalPull,gearGrinding,forwardPull,reversePull},PlotLegends->Placed[PointLegend[{Style["total",Bold,textSize-6,FontFamily->"Helvetica"],Style["gear\ngrinding",Bold,textSize-6,FontFamily->"Helvetica"],Style["reverse cogs\ncontribution",Bold,textSize-6,FontFamily->"Helvetica"],Style["forward cogs\ncontribution",Bold,textSize-6,FontFamily->"Helvetica"]},LegendMarkerSize->textSize-6,LegendMarkers->"Filled"],Right],Joined->True,PlotRange->All,Frame->{{True,False},{True,False}},FrameLabel->{Style["\[Xi] [degrees]",Bold,textSize,FontFamily->"Helvetica"],Style["kcal/mol",textSize,Bold]},PlotStyle->{Darker[Green],Purple,Orange,Blue},FrameTicks->{{Automatic,None},{{-180,-160,-140,-120,-100,-80,-60},None}},LabelStyle->Directive[Black,textSize-6,Bold,FontFamily->"Helvetica"],AspectRatio->0.8,ImageSize->size,TicksStyle->Directive[textSize-6,Bold,FontFamily->"Helvetica"]];
Export[prefix<>"PMF.pdf",p8];
Print[p8];

{p1,p2,p3,p4,p5,p6,p8}
  ]


End[];
EndPackage[];
