(* ::Package:: *)

<<Combinatorica`


Print["Wait...Loading AG..."];


(* Implementacion...
GApops = 
Table[
  Evolution[
    GA[20,           (*Tamanio de la poblacion*)
	  2,             (*Numero de individuos involucrados en la cruza*)
	  PLUS,          (*PLUS: Se seleccionan individuos tanto de los padres como de los hijos para formar las siguiente generacion
					   COMMA: Se descartan los padres y solo los hijos pasan a la siguienete generacion*)
	  10,             (*Numero de hijos*)
	  10,             (*Numero de generaciones*)

      ChromosomeWidth -> 12,
      Partitions -> 2,
      FitnessInterval -> {0,100},
      SelectionMode -> FITPROP,
    
      MutationProbability -> .5,
      PointMutationProbability -> .2,

      RecombinationProbability -> .5,
      RecombinationMode -> POINT[1],
      Homologous -> True,

      InversionProbability -> 0,
      DeletionProbability -> 0,
      DuplicationProbability -> 0,

      EvaluationFunction -> evalFunc,
 
      InitialPopulation -> {}

    ]
  ]//First,
  {3}
];*)


(* ::Subsubsection::Closed:: *)
(*Seleccion Aleatoria de los elementos del alfabeto....*)


randomSelect[s_List] := 
  s[[RandomInteger[{1,Length[s]}]]];


randomSelectMultiple[s_List,n_Integer:1,opts___] :=
	Module[{m,i},
		m=Map[s[[#]]&,i=Take[RandomPermutation[Length[s]],n]];
		If[ReturnIndices /. {opts} /. Options[randomSelectMultiple],
			{m,i},
			m
		]
	]


Options[randomSelectMultiple] = {ReturnIndices -> False};


(* ::Subsubsection::Closed:: *)
(*Generar cromosomas ...*)


createChromosome[width_:1,opts___] :=
	Module[{alph},
	alph := Alphabet /. {opts} /. Options[createChromosome];
	If[RandomLoci /. {opts} /. Options[createChromosome],
		indices = RandomPermutation[width],
		indices = Range[width]];
	chromo[q[undef],p @@ Map[{#,Table[randomSelect[alph],
		{PolyPloidy /. {opts}/. Options[createChromosome]}]}&,indices]]
];


Options[createChromosome] ={
	Alphabet -> {0,1},
	PolyPloidy -> 1,         (*Numero de ramas (= tamanio) del cromosoma*)
	RandomLoci -> False      (*Determina el orden para la clasificacion de los genes
								FALSE: Los genes se ordenaran de acuerdo a su indice
								TRUE: Lo genes se ordenaran aleatoriamente*)
};


(* ::Subsubsection::Closed:: *)
(*Compact Form...*)


ClearAll[CompactForm]
ClearAll[CompactTableForm]


CompactForm[chromo[qual_q,params_p],opts___] :=
  chromo[
    qual,
    If[SortLoci /. {opts} /. Options[CompactForm],
      p @@ Transpose[Map[Last,List @@ Sort[params]]],
      p @@ Transpose[Map[Last,List @@ params]]
    ]
  ]


CompactForm[{chromos__chromo}, opts___] :=
  Map[CompactTableForm[#,opts]&, {chromos}]


Options[CompactForm] = {SortLoci -> False};


CompactTableForm[chromo[qual_q,params_p],opts___] :=
  chromo[
    qual,
    If[SortLoci /. {opts} /. Options[CompactTableForm],
      p[TableForm[Map[Last,List @@ Sort[params]]//Transpose]],
      p[TableForm[Map[Last,List @@ params]//Transpose]]
    ]
  ]


CompactTableForm[{chromos__chromo}, opts___] :=
  Map[CompactTableForm[#,opts]&, {chromos}]


Options[CompactTableForm] = {SortLoci -> False};


(* ::Subsubsection::Closed:: *)
(*Crear una poblacion ...*)


createChromosomes[howMany_Integer:1,opts___] := pop @@ Table[
	createChromosome[ChromosomeWidth /. {opts} /. Options[createChromosomes],opts],
    {howMany}]


Options[createChromosomes] = {
	ChromosomeWidth -> 1, Alphabet -> {0,1}
};


(* ::Subsubsection::Closed:: *)
(*Extraccion de parametros ...*)


(*Del cromosoma... *)
extractParams[chromo[_q,params_p]] := Last /@ params


(* Del vector de parametros... *)
extractParams[params_p] := Last /@ params


(* ::Subsubsection::Closed:: *)
(*Visualizar cromosomas ...*)


(* ::Text:: *)
(*Representacion idividual de cromosomas*)


chromosomePlot[c_chromo,opts___] := Module[
	{params,indices,l,w,alph,trans},
	If[SortLoci /. {opts} /. Options[chromosomePlot],
		params = List @@ Map[Last,Sort[c[[2]]]];
		indices = Range[Length[params]],
		params = List @@ Map[Last,c[[2]]];
		indices = List @@ Map[First,c[[2]]]
	];
	l = Length[params]; (* Number of loci *)
	w = Length[Transpose[params]]; (* chromo strands *)
	alph = Alphabet /. {opts} /. Options[chromosomePlot];
	trans = MapIndexed[#1 -> First[#2]&,alph];
	ListDensityPlot[{Table[-1,{l}], Sequence @@ 
		Reverse[Transpose[params /. trans]],
		Table[-1,{l}]},
		Sequence @@ DeleteCases[{opts},
			(Alphabet -> _) | (SortLoci -> _)],
		Mesh -> True,
		FrameTicks -> {MapIndexed[{First[#2]-.5,#1}&,indices],
		MapIndexed[{First[#2]+.5,#1}&,Range[w,1,-1]]},
		DisplayFunction -> $DisplayFunction
	]
]


Options[chromosomePlot] = {
	Alphabet -> {0,1},
	SortLoci -> False
};


(* ::Text:: *)
(*Visualizacion de miltiples cromosomas...*)


chromosomePictogramPlot[chromos_List, opts___] :=
Module[{},

  MapIndexed[
    chromosomePictogramPlot[#,
      (*PlotLabel -> "Chromosome " <> ToString[First[#2]],*)
      PlotLabel -> None,
      opts]&,
    chromos]
]


(* ::Text:: *)
(*Chromosome Plot*)


chromosomePictogramPlot[chromo[_q, params_p], opts___] :=
Module[{pureParams, l, mp = Length[Last[First[params]]],
        indices},

  indices    = List @@ Map[First,params];
  pureParams = Transpose[List @@ Map[Last,params]];
  l = Length[pureParams];
  
  chromoPlots =
  MapThread[
    pictoPlot[#1, #2, Index -> #3, Sequence @@
      If[GeneIndices /. {opts} /. Options[chromosomePictogramPlot],
        {GeneIndices -> True, Indices -> indices},
        {GeneIndices -> False}
      ],
      opts, Sequence @@ Options[chromosomePictogramPlot]]&,
    {pureParams, 
     Join[
       Table[AxesStyle -> GrayLevel[1],{l-1}],
       {AxesStyle -> GrayLevel[1]}
     ],
     Range[1,l]
    }
  ];
  
  Show[GraphicsArray[List /@ chromoPlots],
    PlotLabel -> (PlotLabel /. {opts} /. {PlotLabel ->
                  FontForm[
                    ToString[mp] <> "-ploid GA-Chromosome, " <>
                    ToString[Length[pureParams//First]] <> 
                    " genes",
                    {"Helvetica",9}]}),
    PlotRange -> All,
    Axes -> False, Frame -> False,
    DisplayFunction -> $DisplayFunction
  ]
]


Options[chromosomePictogramPlot] :=
  Join[
    Options[pictoPlot], 
    Options[pictoMap]
  ]


(* ::Text:: *)
(*Multiple Chromosome Plot*)


chromosomePictogramMultiplePlot[chromos_List, opts___] :=
Module[{pureParams, l, mp, offset, offsetScale, tw, indices,
        graphicsOpts, recMasks, linePlots = {}, chromoPlots,
        lineColFun, lineCol, lineThickn, lineColDec},

  graphicsOpts = GraphicsOptions /. {opts} 
                 /. {GraphicsOptions -> {}};
  
  params = Map[Last,chromos];
  
  mp = Length[Last[First[params]]];
  
  If[Cases[{opts}, TableWidth, Infinity] === {},
    offsetScale = 1; tw = Length[First[params]],
    offsetScale = 
      Ceiling[
        Length[First[params]]/(tw = (TableWidth /. {opts}))]
  ];
  
  chromoPlots =
  MapIndexed[
   (pureParams = Transpose[List @@ Map[Last,#1]];
    l = Length[pureParams];
    offset = 2 - offsetScale First[#2];
    
    MapThread[
      pictoPlot[#1, #2, Index -> #3, 
        VerticalOffset -> offset,
        Sequence @@ DeleteCases[{opts},GeneIndices -> _], 
        Sequence @@ Options[chromosomePictogramPlot]]&,
      {pureParams, 
       Join[
         Table[AxesStyle -> GrayLevel[0],{l-1}],
         {AxesStyle -> GrayLevel[1]}
       ],
       Range[1,l]
      }
    ])&,
    params
  ];
  
  recMasks = RecombinationMasks /. {opts} 
             /. Options[chromosomePictogramMultiplePlot];
             
  If[recMasks != {},
  
    lineColFun = LineColorFunction /. {opts} 
                 /. Options[chromosomePictogramMultiplePlot];
    lineCol    = LineColor /. {opts} 
                 /. Options[chromosomePictogramMultiplePlot];
    lineColDec = LineColorDecrement /. {opts} 
                 /. Options[chromosomePictogramMultiplePlot];
    lineThickn = LineThickness /. {opts} 
                 /. Options[chromosomePictogramMultiplePlot];
    lineCol += lineColDec;
       
    linePlots = Map[
      Show[Graphics[
       {lineColFun[lineCol = Mod[lineCol -lineColDec,1]],
        Thickness[lineThickn],
        Line[
          lin = MapIndexed[
                  getCoord[{#1,Mod[First[#2],tw+1]},tw]&,
                  #
            ]]}],
        DisplayFunction -> Identity
      ]&,
      recMasks
    ]
  ];
  
  Show[If[linePlots === {},
         {chromoPlots},{linePlots,chromoPlots}] // Sequence,
    DefaultFont -> {"Helvetica-Bold",9},
    Frame -> False,
    Axes -> 
    If[GeneIndices /. {opts} 
       /. Options[chromosomePictogramMultiplePlot],
      {True,False}, 
      False
    ],
    AxesOrigin -> {0,-mp},
    AxesStyle -> GrayLevel[0],
    Ticks -> 
    If[GeneIndices /. {opts} 
       /. Options[chromosomePictogramMultiplePlot],
       {List @@ MapIndexed[{First[#2],#1}&,
                 Map[First,params//First]],{}}, 
       None
    ],
    AspectRatio -> 
      (AspectRatio /. graphicsOpts 
       /. Options[chromosomePictogramMultiplePlot]),
    DisplayFunction -> $DisplayFunction
  ]
]


Options[chromosomePictogramMultiplePlot] =
{
  RecombinationMasks -> {},
  AspectRatio -> Automatic,
  LineThickness -> .015,
  LineColorFunction -> Hue,
  LineColor -> .9,
  LineColorDecrement -> .2
};


getCoord[{i1_,i2_},width_] := 
{Mod[i2,width+1],-(i1-1) Ceiling[i2/width]}


(* ::Text:: *)
(*Pictogram Plot on Discrete Lists*)


pictoPlot[params_List,opts___] :=
Module[{width, height, parLength, parList, m, coords, indices,
        hs, vs, textPlot, circlePlot, alph,
        vertOffset, horiOffset},
  
  vertOffset = VerticalOffset /. {opts} /. Options[pictoPlot];
  horiOffset = HorizontalOffset /. {opts} /. Options[pictoPlot];
  
  parLength = Length[params];

  width = TableWidth /. {opts} 
          /. {TableWidth -> parLength};
          
  If[(m = Mod[parLength,width]) == 0,
    parList = params,
    parList = Join[params,Table[" ",{width-m}]];
  ];
 
  alph = Alphabet /. {opts} 
         /. {Alphabet -> DeleteCases[Union[parList]," "]};
         
  indices = 
  If[GeneIndices /. {opts} /. Options[pictoPlot],
    Take[Indices /. {opts} /. Options[pictoPlot],width],
    {}
  ];

  height = Floor[Length[parList]/width];  
  
  coords = 
  Map[{horiOffset,vertOffset} + {1,-1} #&,
    Flatten[Transpose[Array[List,{width,height}]],1]];
  
  textPlot =
  TextListPlot[
    MapThread[Flatten[{#1,#2}]&,{coords,parList}],
    Axes -> False,
    PlotRange -> All,
    DefaultFont -> {"Helvetica-Bold",9},
    DisplayFunction -> Identity
  ];
  
  circlePlot =
  Show[
    Graphics[
      pictoMap[parList,coords,alph,.4,opts]
    ],
    DisplayFunction -> Identity
  ];
  
  minYCoord = Min[Map[Last,coords]];
  
  Show[circlePlot,textPlot,
    GraphicsOptions /. {opts} /. Options[pictoPlot],
    PlotRange -> All,
    DefaultFont -> {"Helvetica-Bold", 9},
    Axes -> 
      If[GeneIndices /. {opts} /. Options[pictoPlot],
        {True, False},
        False
      ], 
    AxesStyle -> 
      {AxesStyle /. {opts} /. {AxesStyle -> GrayLevel[0]},
       Automatic},
    Ticks -> 
      If[GeneIndices /. {opts} /. Options[pictoPlot],
        {MapThread[{#1,#2}&,
          {Range[Length[indices]],indices}],{}},
        None
      ],
    AxesOrigin -> {0,minYCoord-0.5},
    Frame -> True, FrameStyle -> GrayLevel[1],
    FrameTicks -> None,
    FrameLabel -> 
      {None,
       FontForm["-" <> ToString[Index /. {opts}] <> "-",
         {"Courier",9}]}, 
    RotateLabel -> False,
    DisplayFunction -> Identity
  ]
]


Options[pictoPlot] =
{
  VerticalOffset -> 0, HorizontalOffset -> 0,
  GeneIndices -> False, Indices -> {},
  GraphicsOptions -> {DisplayFunction -> Identity}
};


(* ::Text:: *)
(*Pictogram Map*)


pictoMap[params_List, coords_List, 
         alphabet_List, scale_:1, opts___] :=
Module[{p,l,h,hDec,colFunc,pictoElements},

  colFunc = ColorFunction /. {opts} /. Options[pictoMap]; 
  h = StartColor /. {opts} /. Options[pictoMap]; 
  hDec = ColorDecrement /. {opts} /. Options[pictoMap];
  p = pictoElements =
    Pictograms /. {opts} /. Options[pictoMap];
    
  l = Length[p];
  
  pictoRules = 
  MapIndexed[
    (If[Mod[First[#2],l] == 0, h = Mod[h - hDec,1]];
    #1 -> {h,First[p = RotateLeft[p]]})&,
    alphabet
  ] ~Join~ {" " -> {1,NULLDISK}};
  
  MapThread[#2[[2]][#1,scale,#2[[1]],colFunc]&,
    {coords,params /. pictoRules}
  ] // Flatten
]


Options[pictoMap] =
{
  Pictograms -> {TRIANGLE, DISK, RECTANGLE, DIAMOND},
  ColorFunction -> Hue,
  StartColor -> 0.9,
  ColorDecrement -> .1
};


DISK[{x_,y_}, s_:.5, g_:0, c_:Hue] := {c[g], Disk[{x,y},s]}

RECTANGLE[{x_,y_}, s_:.5, g_:0, c_:Hue] := 
  {c[g], Rectangle[{x-s,y-s},{x+s,y+s}]}
                              
DIAMOND[{x_,y_}, s_:.5, g_:0, c_:Hue] := 
  {c[g], Polygon[{{x,y-s},{x-s,y},{x,y+s},{x+s,y}}]}
  
TRIANGLE[{x_,y_}, s_:.5, g_:0, c_:Hue] := 
  {c[g], Polygon[{{x-s,y-s},{x,y+s},{x+s,y-s}}]}
  
NULLDISK[{x_,y_}, s_, g_, c_] := DISK[{x,y},0]


(* ::Text:: *)
(*Faces*)


(* ::Text:: *)
(*Visualizing Faces*)


Off[General::spell1]


face::usage = "face[{headecc_, eyesize_, eyespacing_, 
		eyeeccent_, pupsize_, browslant_, nozesize_,
		mouthshape_, mouthsize_, mouthopening_}]";


face[{headecc_, eyesize_, eyespacing_, eyeeccent_, 
	  pupsize_, browslant_, nozesize_, mouthshape_,
	  mouthsize_, mouthopening_}] :=
    Graphics[Flatten[
       {head[headecc],
        eyes[eyesize, eyespacing, eyeeccent, 
        	 2 pupsize],
        brows[browslant], nose[nozesize],
        mouth[mouthshape, mouthsize, mouthopening]}], 
                    AspectRatio -> Automatic,  
                    PlotRange-> {{-1.5, 1.5}, 
                    			 {-1.5, 1.5}}]


head[eccent_] := Block[
 		{xrad = 1 + (eccent - 5) / 25, 
  		 yrad = 1 - (eccent - 5) / 25},
  		 {GrayLevel[.9],
          Disk[{0, 0}, (1 + Abs[eccent - 5]/20)
            {xrad, yrad}],
          GrayLevel[0],
          Circle[{0, 0}, (1 + Abs[eccent - 5]/20) 
            {xrad, yrad}]}]


eyes[size_, spacing_, eccent_, pupsize_] := Block[
        {xcenter = (1/3) + (spacing - 5) / 30,
         xrad = (1/6) + ((size - 5) + (eccent - 5)) 
         				  / 70,
         yrad = (1/6) + ((size - 5) - (eccent - 5)) 
         				  / 70},
  {GrayLevel[.8],
    Disk[{xcenter, 1/3},  {xrad, yrad}],
   GrayLevel[0],
    Circle[{xcenter, 1/3},  {xrad, yrad}],
    PointSize[(pupsize + 1) / 200], 
    Point[{xcenter, 1/3}],
   GrayLevel[.8],
    Disk[{-xcenter, 1/3},  {xrad, yrad}],
   GrayLevel[0],
    Circle[{-xcenter, 1/3},  {xrad, yrad}],
    PointSize[(pupsize + 1) / 200], 
    Point[{-xcenter, 1/3}]
  }
                                                  ]


brows[slant_] :=  Block[
       {xstart = (1/3) - (1/6) Cos[(slant - 5) Pi 
       							    / 20],
        ystart = (2/3) - (1/6) Sin[(slant - 5) Pi 
        							/ 20],
        xend = (1/3) + (1/6) Cos[(slant - 5) Pi 
        						  / 20],
        yend =  (2/3) + (1/6) Sin[(slant - 5) Pi 
        						   / 20]},
    {Line[{{xstart, ystart}, {xend, yend}}],
     Line[{{-xstart, ystart}, {-xend, yend}}]}
                       ]


nose[size_] := Block[
        {scale = 1 + (size - 5) / 13},
 {
   GrayLevel[0.6],
   Polygon[scale {{0, 1/6}, {-1/6, -1/6}, 
              {1/6, -1/6}, {0, 1/6}}],
   GrayLevel[0],
   Line[scale {{0, 1/6}, {-1/6, -1/6}, 
              {1/6, -1/6}, {0, 1/6}}]
 }
                    ]


mouth[shape_, size_, opening_] := Block[
  {fx, gx, xstart, xend, ystart, ymin, ymax, xstep },
 xstart = -1/3 - (size - 5) / 15; 
 xend = 1/3 + (size - 5) / 15;
 ystart = -1/2 + (shape - 5) * size / 150;
 ymax = -1/2 + (.9 opening - 1) / 27;
 ymin = -1/2 - (.9 opening - 1) / 30;
 fx = Fit[{{xstart, ystart}, {0, ymax}, 
 		   {xend, ystart}}, {1, x, x^2}, x];
 gx = Fit[{{xstart, ystart}, {0, ymin},
 		   {xend, ystart}}, {1, x, x^2}, x];
 xstep = (xend - xstart) / 10;
    {GrayLevel[1],
     Polygon[Table[{x, fx}, {x, xstart, xend, xstep}]],
     Polygon[Table[{x, gx}, {x, xstart, xend, xstep}]],
     GrayLevel[0],
     Line[Table[{x, fx}, {x, xstart, xend, xstep}]],
     Line[Table[{x, gx}, {x, xstart, xend, xstep}]]}
                                        ]


On[General::spell1]


(* ::Text:: *)
(*GA Chromosome Representation by Faces*)


chromosomeFacePlot[chromo[_q,params_p], opts___] :=
Module[{},

  parList = Flatten[List @@ Map[Last,params]];
  
  Show[face[parList], 
    Sequence @@ DeleteCases[{opts}, TableWidth -> _]
  ]
]


chromosomeFacePlot[chromos_List, opts___] :=
Module[{tw, faceGraphics},

  faceGraphics =
  Map[
    chromosomeFacePlot[#,DisplayFunction -> Identity,opts]&,
    chromos];
    
  If[Cases[{opts},TableWidth,2] === {},
  
    Show[GraphicsArray[faceGraphics],
      DisplayFunction -> $DisplayFunction
    ],
    
    Show[GraphicsArray[
      Partition[faceGraphics,
        TableWidth /. {opts} /. Options[chromosomeFacePlot]
      ]],
      DisplayFunction -> $DisplayFunction
    ]
  ];

]


(* ::Subsubsection::Closed:: *)
(*GA Chromosomes Decoding*)


Decoding[chromo[qual_q,params_p],opts___] :=
Module[{paramsList},
  
  paramsList = 
  Partition[Flatten[List @@ extractParams[params]],
    Floor[Length[params]/(Partitions /. {opts} 
      /. Options[Decoding])]];
  
  chromo[
    qual,
    p @@ Map[DecodingFunction /. {opts} 
      /. Options[Decoding],
      If[(Encoding /. {opts} /. Options[Decoding])
         === GRAY,
         grayToStandard /@ paramsList,
         paramsList
      ]
    ]
  ]
]


Decoding[params_p,opts___] :=
Module[{paramsList},

  paramsList = 
  Partition[Flatten[List @@ extractParams[params]],
    Floor[Length[params]/(Partitions /. {opts} 
      /. Options[Decoding])]];
  
  p @@ Map[DecodingFunction /. {opts} /. Options[Decoding],
         If[(Encoding /. {opts} /. Options[Decoding])
            === GRAY,
           grayToStandard /@ paramsList,
           paramsList
         ]
       ]
]


Decoding[population_pop, opts___] :=
Module[{},

  Map[Decoding[#,opts]&, population]
]


Decoding[population_List, opts___] :=
Module[{},

  Map[Decoding[#,opts]&, population]
]


Decoding[x_,___] := x


Options[Decoding] =
{
  Alphabet -> {0,1},
  Encoding -> STANDARD, (* GRAY *)
  Partitions -> 2,
  DecodingFunction -> (lettersToNumber[#,2]&)
};


(* ::Subsubsection:: *)
(*Evaluacion ...*)


Evaluation[chromo[qualities_q,params_p,x___],opts___] := 
	Module[{g},
		g = EvaluationFunction /. {opts} /. Options[Evaluation];
		chromo[q[g[extractParams[params]]//N],params,x]
	]


Options[Evaluation] = {
	EvaluationFunction :> (Apply[Plus,#]&)
};


(* ::Subsubsection::Closed:: *)
(*Mutacion de Cromosomas ...*)


(* ::Text:: *)
(*Utility Function*)


chi[] := RandomReal[{0,1}];


(* ::Text:: *)
(*Mutacion Simple ...*)


Mutation[chromo[qual_q,params_p],opts___] := Module[{alph,flipProb},
	If[chi[] <= (MutationProbability /. {opts} /. Options[Mutation]),
		alph := Alphabet /. {opts} /. Options[Mutation];
		flipProb = PointMutationProbability /. {opts} /. Options[Mutation];
		chromo[q[undef],Map[geneMutation[#,flipProb,alph]&,params]],
		chromo[qual,params]
	]
]


Options[Mutation] ={
	Alphabet -> {0,1},
	MutationProbability -> 1,
	PointMutationProbability -> .5};


(* ::Text:: *)
(*Mutacion Poliploide ...*)


geneMutation[{i_,allels_List},flipProb_:.5,alphabet_:{0,1}] := Module[{},
	{i,Map[If[chi[]<=flipProb,flip[#,alphabet],#]&,allels]}
]


flip[e_,sel_List] := randomSelect[Complement[sel,{e}]]


(* ::Subsubsection::Closed:: *)
(*Recombinacion ...*)


Recombination[chromos_List,opts___] := Module[
	{params, paramsSorted, indices,recParams, recMasks},
	If[chi[] > (RecombinationProbability /. {opts} /. Options[Recombination]),
		Return[{chromos//First,{}}]];               
  (* Extract chromo parameters *)
	params = Map[List @@ Last[#]&,chromos];
  (* Homologous? *)
  If[Homologous /. {opts} /. Options[Recombination],
    recombineHomologous[params, opts],
    recombineNonHomologous[params, opts] (*Aunque esta no se ocupa
										 xq los cromosomas son homologos*)
  ]
]


Options[Recombination] := Join[{RecombinationProbability -> 1,
	Homologous -> True},
	Options[Recombine]
]


(* ::Text:: *)
(*Recombinacion de cromosomas homologos ...*)


recombineHomologous[params_, opts___] := Module[
	{paramsSorted, recParams, recMasks},
	(* Extract indices *)
	indices = Map[Map[First,#]&,params];
	(* Sort parameters according to the first chromo *)
	paramsSorted = Map[Sort[#][[First[indices]]]&,params]; (*ok*)
	(* Get rid of indices *)
	paramsSorted = Transpose[
    Map[Transpose[Map[Last,#]]&,paramsSorted]];
	(* Recombine parameters *)
	{recParams,recMasks} = 
		Transpose[Map[Recombine[#,opts]&,paramsSorted]];
	{chromo[q[undef],p @@ 
		MapThread[{#1,#2}&,{indices//First,recParams//Transpose}]],
		recMasks}
]


(* ::Text:: *)
(*Recombinacion de cromosomas no homologos ...*)


recombineNonHomologous[params_, opts___] := 
	Module[{paramsSorted, recParams, recMasks},
		(* Sort parameters *)
		paramsSorted = Sort /@ params;
		(* Get rid of indices *)
		paramsSorted = Transpose[Map[
			Transpose[Map[Last,#]]&,paramsSorted]];
		(* Recombine parameters *)
		{recParams,recMasks} = Transpose[Map
			[Recombine[#,opts]&,paramsSorted]];
		{chromo[q[undef],p @@
			MapIndexed[{First[#2],{#1}}&,recParams//First]],
			recMasks}
	]


Recombine[params_List,opts___] := Module[
	{rm,coPoints,l = Length[params//First],paramsLength,
	recFun},
		paramsLength = Length[params]; (*ok*)
		Switch[rm = (RecombinationMode /. {opts} /. Options[Recombine]),
			POINT[_],
				coPoints = First[rm];
				{MapIndexed[params[[#1,First[#2]]]&,
					mask = recombinationMask[l,
						RecombinationPoints -> coPoints,
						RecombinationPartners -> Length[params]]
					],
				mask},
			MULTIPOINT,
				{MapIndexed[params[[#1,First[#2]]]&,
					mask = Table[RandomInteger[{1,Length[params]}],{l}]],
				mask},
			MASK[_],
				mask = First[rm];
				{MapIndexed[params[[#1,First[#2]]]&,mask],mask},
			CROSSING[_],
				Crossing[params, Crossings -> First[rm], opts],
			NORMAL|_,
				recFun := RecombinationFunction /. {opts} 
							/. Options[Recombine];
				{MapIndexed[recFun[#1,First[#2]]&,Transpose[params]],{}}
		]
	]


Options[Recombine] := {
	RecombinationMode -> CROSSING[1],
	RecombinationFunction :> r 
};


nonHomologous = {{a,b,c,d,e,f,g,h,k,l},{A,B,C,D,E,F,G,H}};


Recombine[nonHomologous, ReturnPair -> False,
  RecombinationMode -> CROSSING[1]]


(* ::Text:: *)
(*Recombinacion Discreta ...*)


DISCRETE = (randomSelect[#]&);


(* ::Text:: *)
(*Mascara de Recombinacion ...*)


recombinationMask[width_Integer, opts___] := Module[{s,ss,r,t,i},
	coPoints = RecombinationPoints /. {opts} /. 
		Options[recombinationMask];
	rePartners = RecombinationPartners /. {opts} /.
		Options[recombinationMask];
	reRange = RecombinationRange /. {opts} /.
		Options[recombinationMask];
	(* Generate random crossover points *)
	cps = Sort[randomSelectMultiple[Range[1,width-1],coPoints]];
	cps = Join[{0},cps,{width}];
	(* Generate combination mask with flipped indices *)
	s = Range[1,rePartners];
	r = ss = randomSelectMultiple[s,reRange];
	mask = Map[(i = randomSelect[ss];
		ss = Complement[r,{i}];
		Table[i,{#}])&,
		Map[#[[2]]-#[[1]]&,Partition[cps,2,1]] ] // Flatten;
	mask
]


Options[recombinationMask] =
{
  RecombinationRange -> 2,
  RecombinationPartners -> 2,
  RecombinationPoints -> 1
};


inverseRecombinationMask[mask_List] :=
Module[{alph},

  alph = Union[mask];
  
  mask /. MapThread[#1 -> #2 &,{alph,RotateLeft[alph]}]
]


(* ::Text:: *)
(*Recombnacion de cromosomas no homologos con diferentes puntos de cruza*)


Crossing[params_List, opts___] :=

Module[{cros, l = Length[params//First],coPoints,
        coRanges, coRangesOne},
  
  cros = Crossings /. {opts} /. Options[Crossing];
  
  parLengths = Map[Length,params];
  
  coPoints =
  Map[Sort[Table[Random[Integer,{2,#}],{cros}]]&,
    parLengths];
  
  coRanges = MapThread[Partition[Join[{1},#1,{#2+1}],2,1]&,
               {coPoints,parLengths}];
  
  coRangesOne =
  Flatten[
    MapThread[{#1,#2}&,
      {Partition[coRanges[[1]],1,2],
       Partition[RotateLeft[coRanges[[2]]],1,2]}],2];
       
  coRangesTwo =
  Flatten[
    MapThread[{#1,#2}&,
      {Partition[coRanges[[2]],1,2],
       Partition[RotateLeft[coRanges[[1]]],1,2]}],2];
  
  If[EvenQ[cros], 
    coRangesOne = Drop[coRangesOne,-1];
    coRangesTwo = Drop[coRangesTwo,-1]
  ];
  
  {If[ReturnPair /. {opts} /. Options[Crossing],
    {Flatten[
     MapIndexed[Take[params[[2-Mod[First[#2],2]]],#1-{0,1}]&,
       coRangesOne],1],
     Flatten[
     MapIndexed[Take[params[[1+Mod[First[#2],2]]],#1-{0,1}]&,
       coRangesTwo],1]},
    Flatten[
     MapIndexed[Take[params[[2-Mod[First[#2],2]]],#1-{0,1}]&,
       coRangesOne],1]
   ],
   coPoints}
]


Options[Crossing] =
{
  Crossings -> 1,
  ReturnPair -> False
};


(* ::Text:: *)
(*Recombinacion Meiotica de cromosomas diploides...*)


meiosis[chromos_genome, opts___] :=
Module[{haploidStrand, params, hapChromoOne, hapChromoTwo,
        recChromoOne, recChromoTwo, recMask},

  haploidStrand[p_,i_] := Map[Map[{#[[1]],{#[[2,i]]}}&,#]&,p];
  
  (* Extract haploid chromosome strands *)
  
  params = Map[Last,chromos];
  
  {hapChromoOne,hapChromoTwo} = 
    Map[chromo[q[undef],haploidStrand[params,#]//First]&,
    {1,2}];
  
  (* Recombine strands via one-point crossover *)
  
  {recChromoOne,recMask} =
    Recombination[{hapChromoOne,hapChromoTwo}, 
      RecombinationMode -> POINT[1]];
    
  (* Recombine with inverse mask *)
  
  {recChromoTwo,recMask} =
    Recombination[{hapChromoOne,hapChromoTwo}, 
      RecombinationMode -> 
        MASK[inverseRecombinationMask[recMask//First]]];
  
  (* Form two genomes with one recombined strand 
     and one original strand *)
     
  {genome[{mergeChromosomes[{hapChromoOne,recChromoOne}]}],
   genome[{mergeChromosomes[{hapChromoTwo,recChromoTwo}]}]}
]


(* ::Text:: *)
(*Separating Polyploid Chromosomes into Haploid Single Strands*)


splitChromosome[chr_chromo, opts___] :=
Module[{haploidStrand, params, ploidy, hapChromos},
  
  haploidStrand[p_,i_] := Map[{#[[1]],{#[[2,i]]}}&,p];
  
  params = chr//Last;
  
  ploidy = Length[Last[params//First]];
  
  hapChromos = 
    Map[chromo[q[undef],haploidStrand[params,#]]&,
      Range[ploidy]]
]


(* ::Text:: *)
(*Mergin Several Polyploid GA Chromosomes into a Single Polyploid GA Chromosome*)


mergeChromosomes[chromos_List, opts___] :=
Module[{params, indices, paramsSorted},

  (* Sort each chromosome genes according to first chromo *)
  
  params = Map[List @@ Last[#]&,chromos];
  
  indices = Map[Map[First,#]&,params];
  
  (* Sort parameters according to the first chromo *)
  
  paramsSorted = Map[Sort[#][[First[indices]]]&,params];
  paramsSorted = Map[Last /@ #&,paramsSorted];
  
  (* Merge genes *)
  
  chromo[q[undef],
    p @@ MapThread[{#2,#1}&,
      {Flatten /@ Transpose[paramsSorted],indices//First}]
  ]
]


(* ::Subsubsection::Closed:: *)
(*Inversion*)


Inversion[chromo[qual_q,params_p], opts___] :=
	Module[{l = Length[params], i1, i2, perm},
		If[chi[] > (InversionProbability /. {opts} /. Options[Inversion]),
			Return[{chromo[qual,params],{}}]
		];
		perm = InvPermutation /. {opts} /. Options[Inversion];
		{i1,i2} = Sort[Table[Random[Integer,{1,l}],{2}]];
		{chromo[q[undef],
			Join[
				Take[params, Max[0,i1-1]],
				perm[Take[params, {i1,i2}]],
				Take[params, Min[l,i2]-l]
			]
		],{i1,i2}}
	]


Options[Inversion] = {
	InvPermutation -> Reverse,
	InversionProbability -> 1
	};


(* ::Subsubsection::Closed:: *)
(*Borrado*)


Deletion[chromo[qual_q,params_p], opts___] :=
	Module[{l = Length[params], low, high},
		If[(chi[] > (DeletionProbability /. {opts} /. Options[Deletion])) ||
			l <= (MinChromoWidth /. {opts} /. Options[Deletion]),
			Return[{chromo[qual,params],{}}]
		];
		{low,high} = Sort[Table[Random[Integer,{1,l}],{2}]];
		If[l-(high-low) <= (MinChromoWidth /. {opts} /. Options[Deletion]),
			Return[{chromo[qual,params],{}}]
		];
		{chromo[q[undef],
			Join[Take[params, Max[0,low-1]],
			Take[params, Min[l,high]-l]]
			],{low,high}}
		]


Options[Deletion] = {
	DeletionProbability -> 1,
	MinChromoWidth -> 2
	};


(* ::Subsubsection::Closed:: *)
(*Duplicado*)


Duplication[chromo[qual_q,params_p], opts___] := 
	Module[{l = Length[params], low, high},
		If[chi[] > (DuplicationProbability /. {opts} 
			/. Options[Duplication]),
			Return[{chromo[qual,params],{}}]
		];
		{low,high} = Sort[Table[Random[Integer,{1,l}],{2}]];
		{chromo[q[undef],
			Join[
				Take[params, Max[0,low-1]],
				Take[params, {low,high}],
				Take[params, {low,high}],
				Take[params, Min[l,high]-l]
				]
			],{low,high}}
		]


Options[Duplication] = {
	DuplicationProbability -> 1
	};


(* ::Subsubsection::Closed:: *)
(*Seleccion*)


(* ::Text:: *)
(*Seleccion en listas*)


ClearAll[selectFitProp]


selectFitProp[fitnesses_, opts___] := Module[{cumFits, fits, r, m},
	fits = If[Negative[m = Min[fitnesses]], fitnesses -m, fitnesses];
	cumFits = FoldList[Plus, 0, fits] // Rest;
	r = Random[] Last[cumFits];
	Position[cumFits, _?(# >= r &), {1}, 1][[1,1]]
]


(* ::Text:: *)
(*Seleccion basada en el rango*)


ClearAll[selectRankBased]


selectRankBased[fitnesses_, opts___] := Module[{fitsAndIndsSorted, index},
	fitsAndIndsSorted = Sort[MapIndexed[{#1,First[#2]}&,fitnesses]];
	index = selectFitProp[Range[Length[fitnesses]]];
	fitsAndIndsSorted[[index]]//Last
]


(* ::Text:: *)
(*Seleccion Elitista*)


selectElite[fitnesses_, opts___] := Module[{bestCount, fitsAndIndsSorted},
	bestCount = Elite /. {opts} /. Options[selectElite];
	fitsAndIndsSorted = 
		Reverse[Sort[MapIndexed[{#1,First[#2]}&,fitnesses]]];
	selectRandom[Last /@ 
		If[IntegerQ[bestCount],
			Take[fitsAndIndsSorted, bestCount],
			Take[fitsAndIndsSorted, 
				Max[1, Floor[bestCount Length[fitnesses]]]]
        ]
     ]
]


Options[selectElite] ={
	Elite -> 1
};


(* ::Text:: *)
(*Seleccion del mejor ...*)


selectBest[fitnesses_, opts___] := Module[{bestCount, fitsAndIndsSorted},
	bestCount = Best /. {opts} /. Options[selectBest];
	fitsAndIndsSorted = Reverse[Sort[MapIndexed[{#1,First[#2]}&,
		fitnesses]]];
	Last /@ If[IntegerQ[bestCount],Take[fitsAndIndsSorted, bestCount],
		Take[fitsAndIndsSorted,
			Max[1, Floor[bestCount Length[fitnesses]]]]
	]
]


Options[selectBest] = {
	Best -> 1
};


(* ::Text:: *)
(*Funcion de seleccion ...*)


GASelection[population_pop,howMany_Integer:1,opts___] := Module[{},
	Switch[SelectionMode /. {opts} /. Options[GASelection],
		BEST,
			Take[population // Sort // Reverse, howMany],
		RANDOM,
			pop @@
			Table[population[[RandomInteger[{1,Length[population]}]]],
			{howMany}],
		FITPROP,
			Table[population[[
				selectFitProp[extractQuals[population],opts]]],
				{howMany}],
		RANKBASED,
			Table[population[[
				selectRankBased[extractQuals[population],opts]]],
				{howMany}],
		ELITE,
			Table[population[[
				selectElite[extractQuals[population],opts,
					Elite -> 1.]]],
				{howMany}]
	]
]


Options[GASelection] = {
	SelectionMode -> BEST (* RANDOM or FITPROP *)
};


(* ::Text:: *)
(*Extraccion de valores de fitness ...*)


extractQuals[population_pop] := List @@ Map[First[First[#]]&,population]


(* ::Subsubsection::Closed:: *)
(*Esquema Evolutivo*)


Evolution[
	GA[m_?IntegerQ,    (*Tamanio de la poblacion*)
	r_?IntegerQ,       (*Tamanio matrimonios*)
	s_?AtomQ,          (*Tipo de seleccion*)
	l_?IntegerQ,       (*# de hijos*)
	g_?IntegerQ,       (*# de Generaciones Permitidas*)
	opts___] i_:1]:=
	Module[{initialParents,evalFct,parents,children,selPool},
	Print["GA evolution ..."];

	(* Generate m initial individuals *)
		If[(initialParents=(InitialPopulation/.{opts}/.Options[Evolution]))==={},
			initialParents=createChromosomes[m,ChromosomeWidth->(ChromosomeWidth/.{opts}/.Options[Evolution]),opts]
			];

			(* Get the evaluation function *)
			evalFct=EvaluationFunction/.{opts}/.Options[Evaluation];

			(* Evaluate initial individuals *)
			initialParents=Map[Evaluation[#,EvaluationFunction:>evalFct]&,initialParents];

			(* Iterate for i independent runs *)
			Table[
				(*Iterate for g generations*)
				NestList[(parents=#;
					children=Apply[pop,ComposeList[{Map[First[Recombination[#,opts]]&,#]&,Map[Mutation[#,opts]&,#]&,
						Map[First[Inversion[#,opts]]&,#]&,Map[First[Deletion[#,opts]]&,#]&,Map[First[Duplication[#,opts]]&,#]&,
						Map[Evaluation[#,EvaluationFunction:>evalFct,opts]&,#]&},
						Partition[GASelection[parents,r*l,SelectionMode->(SelectionMode/.{opts}/.Options[Evolution]),opts],r]
						]//Last
					];
					selPool = Switch[s,
						PLUS,Join[parents,children],
						COMMA,children];
					GASelection[selPool,m,SelectionMode->BEST])&,initialParents,g
				],(*end nesting*)
				{i} (*independent runs*)
			]
		]


Options[Evolution] = 
	Join[{
		InitialPopulation -> {},
		 SelectionMode -> FITPROP},
		Options[createChromosomes],
		Options[Evaluation],
		Options[Mutation],
		Options[Recombination],
		Options[Inversion],
		Options[Deletion],
		Options[Duplication]
	];
