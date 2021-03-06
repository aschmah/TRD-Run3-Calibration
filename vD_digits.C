void vD_digits()
{
//=========Macro generated from canvas: tg_vdrift_vs_det_can/tg_vdrift_vs_det_can
//=========  (Tue Apr  7 16:07:29 2020) by ROOT version 6.18/04
   TCanvas *tg_vdrift_vs_det_can = new TCanvas("tg_vdrift_vs_det_can", "tg_vdrift_vs_det_can",59,121,900,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   tg_vdrift_vs_det_can->Range(-144,-0.8571429,576,3.428571);
   tg_vdrift_vs_det_can->SetFillColor(10);
   tg_vdrift_vs_det_can->SetBorderMode(0);
   tg_vdrift_vs_det_can->SetBorderSize(2);
   tg_vdrift_vs_det_can->SetTickx(1);
   tg_vdrift_vs_det_can->SetTicky(1);
   tg_vdrift_vs_det_can->SetLeftMargin(0.2);
   tg_vdrift_vs_det_can->SetRightMargin(0.05);
   tg_vdrift_vs_det_can->SetBottomMargin(0.2);
   tg_vdrift_vs_det_can->SetFrameFillColor(0);
   tg_vdrift_vs_det_can->SetFrameBorderMode(0);
   tg_vdrift_vs_det_can->SetFrameBorderMode(0);
   
   TH1D *tg_vdrift_vs_det_candummy_copy__1 = new TH1D("tg_vdrift_vs_det_candummy_copy__1","",1000,0,540);
   tg_vdrift_vs_det_candummy_copy__1->SetMinimum(0);
   tg_vdrift_vs_det_candummy_copy__1->SetMaximum(3);
   tg_vdrift_vs_det_candummy_copy__1->SetDirectory(0);
   tg_vdrift_vs_det_candummy_copy__1->SetStats(0);
   tg_vdrift_vs_det_candummy_copy__1->SetLineColor(10);
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->SetTitle("det");
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->SetRange(1,1000);
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->CenterTitle(true);
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->SetNdivisions(505);
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->SetLabelFont(42);
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->SetLabelSize(0.06);
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->SetTitleSize(0.06);
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->SetTitleOffset(1.2);
   tg_vdrift_vs_det_candummy_copy__1->GetXaxis()->SetTitleFont(42);
   tg_vdrift_vs_det_candummy_copy__1->GetYaxis()->SetTitle("drift velocity");
   tg_vdrift_vs_det_candummy_copy__1->GetYaxis()->CenterTitle(true);
   tg_vdrift_vs_det_candummy_copy__1->GetYaxis()->SetNdivisions(505);
   tg_vdrift_vs_det_candummy_copy__1->GetYaxis()->SetLabelFont(42);
   tg_vdrift_vs_det_candummy_copy__1->GetYaxis()->SetLabelSize(0.06);
   tg_vdrift_vs_det_candummy_copy__1->GetYaxis()->SetTitleSize(0.06);
   tg_vdrift_vs_det_candummy_copy__1->GetYaxis()->SetTitleOffset(1.2);
   tg_vdrift_vs_det_candummy_copy__1->GetYaxis()->SetTitleFont(42);
   tg_vdrift_vs_det_candummy_copy__1->GetZaxis()->SetLabelFont(42);
   tg_vdrift_vs_det_candummy_copy__1->GetZaxis()->SetLabelSize(0.035);
   tg_vdrift_vs_det_candummy_copy__1->GetZaxis()->SetTitleSize(0.035);
   tg_vdrift_vs_det_candummy_copy__1->GetZaxis()->SetTitleOffset(1);
   tg_vdrift_vs_det_candummy_copy__1->GetZaxis()->SetTitleFont(42);
   tg_vdrift_vs_det_candummy_copy__1->Draw("h");
   
   Double_t tg_vdrift_vs_det_fx1[540] = {
   0,
   1,
   2,
   3,
   4,
   5,
   6,
   7,
   8,
   9,
   10,
   11,
   12,
   13,
   14,
   15,
   16,
   17,
   18,
   19,
   20,
   21,
   22,
   23,
   24,
   25,
   26,
   27,
   28,
   29,
   30,
   31,
   32,
   33,
   34,
   35,
   36,
   37,
   38,
   39,
   40,
   41,
   42,
   43,
   44,
   45,
   46,
   47,
   48,
   49,
   50,
   51,
   52,
   53,
   54,
   55,
   56,
   57,
   58,
   59,
   60,
   61,
   62,
   63,
   64,
   65,
   66,
   67,
   68,
   69,
   70,
   71,
   72,
   73,
   74,
   75,
   76,
   77,
   78,
   79,
   80,
   81,
   82,
   83,
   84,
   85,
   86,
   87,
   88,
   89,
   90,
   91,
   92,
   93,
   94,
   95,
   96,
   97,
   98,
   99,
   100,
   101,
   102,
   103,
   104,
   105,
   106,
   107,
   108,
   109,
   110,
   111,
   112,
   113,
   114,
   115,
   116,
   117,
   118,
   119,
   120,
   121,
   122,
   123,
   124,
   125,
   126,
   127,
   128,
   129,
   130,
   131,
   132,
   133,
   134,
   135,
   136,
   137,
   138,
   139,
   140,
   141,
   142,
   143,
   144,
   145,
   146,
   147,
   148,
   149,
   150,
   151,
   152,
   153,
   154,
   155,
   156,
   157,
   158,
   159,
   160,
   161,
   162,
   163,
   164,
   165,
   166,
   167,
   168,
   169,
   170,
   171,
   172,
   173,
   174,
   175,
   176,
   177,
   178,
   179,
   180,
   181,
   182,
   183,
   184,
   185,
   186,
   187,
   188,
   189,
   190,
   191,
   192,
   193,
   194,
   195,
   196,
   197,
   198,
   199,
   200,
   201,
   202,
   203,
   204,
   205,
   206,
   207,
   208,
   209,
   210,
   211,
   212,
   213,
   214,
   215,
   216,
   217,
   218,
   219,
   220,
   221,
   222,
   223,
   224,
   225,
   226,
   227,
   228,
   229,
   230,
   231,
   232,
   233,
   234,
   235,
   236,
   237,
   238,
   239,
   240,
   241,
   242,
   243,
   244,
   245,
   246,
   247,
   248,
   249,
   250,
   251,
   252,
   253,
   254,
   255,
   256,
   257,
   258,
   259,
   260,
   261,
   262,
   263,
   264,
   265,
   266,
   267,
   268,
   269,
   270,
   271,
   272,
   273,
   274,
   275,
   276,
   277,
   278,
   279,
   280,
   281,
   282,
   283,
   284,
   285,
   286,
   287,
   288,
   289,
   290,
   291,
   292,
   293,
   294,
   295,
   296,
   297,
   298,
   299,
   300,
   301,
   302,
   303,
   304,
   305,
   306,
   307,
   308,
   309,
   310,
   311,
   312,
   313,
   314,
   315,
   316,
   317,
   318,
   319,
   320,
   321,
   322,
   323,
   324,
   325,
   326,
   327,
   328,
   329,
   330,
   331,
   332,
   333,
   334,
   335,
   336,
   337,
   338,
   339,
   340,
   341,
   342,
   343,
   344,
   345,
   346,
   347,
   348,
   349,
   350,
   351,
   352,
   353,
   354,
   355,
   356,
   357,
   358,
   359,
   360,
   361,
   362,
   363,
   364,
   365,
   366,
   367,
   368,
   369,
   370,
   371,
   372,
   373,
   374,
   375,
   376,
   377,
   378,
   379,
   380,
   381,
   382,
   383,
   384,
   385,
   386,
   387,
   388,
   389,
   390,
   391,
   392,
   393,
   394,
   395,
   396,
   397,
   398,
   399,
   400,
   401,
   0,
   0,
   0,
   0,
   0,
   0,
   408,
   409,
   410,
   411,
   412,
   413,
   414,
   415,
   416,
   417,
   418,
   419,
   420,
   421,
   422,
   423,
   424,
   425,
   426,
   427,
   428,
   429,
   430,
   431,
   0,
   0,
   0,
   0,
   0,
   0,
   438,
   439,
   440,
   441,
   442,
   443,
   444,
   445,
   446,
   447,
   448,
   449,
   450,
   451,
   452,
   453,
   454,
   455,
   456,
   457,
   458,
   459,
   460,
   461,
   0,
   0,
   0,
   0,
   0,
   0,
   468,
   469,
   470,
   471,
   472,
   473,
   474,
   475,
   476,
   477,
   478,
   479,
   480,
   481,
   482,
   483,
   484,
   485,
   486,
   487,
   488,
   489,
   490,
   491,
   492,
   493,
   494,
   495,
   496,
   497,
   498,
   499,
   500,
   501,
   502,
   503,
   504,
   505,
   506,
   507,
   508,
   509,
   510,
   511,
   512,
   513,
   514,
   515,
   516,
   517,
   518,
   519,
   520,
   521,
   522,
   523,
   524,
   525,
   526,
   527,
   528,
   529,
   530,
   531,
   532,
   533,
   534,
   535,
   536,
   537,
   538,
   539};
   Double_t tg_vdrift_vs_det_fy1[540] = {
   1.499244,
   1.475175,
   1.335216,
   1.473563,
   1.493113,
   0.4704081,
   1.488364,
   1.471577,
   1.148198,
   1.461033,
   1.486204,
   1.470804,
   1.133727,
   1.471219,
   1.457179,
   1.205273,
   1.471699,
   1.335216,
   1.4855,
   1.470656,
   1.470235,
   1.471008,
   1.375497,
   1.478952,
   1.487009,
   1.341701,
   0.4719917,
   1.335216,
   1.359682,
   0.9618259,
   0.4281277,
   1.409255,
   1.409255,
   1.480425,
   1.506014,
   1.479533,
   1.409255,
   1.48723,
   1.476918,
   1.36471,
   1.409255,
   1.480073,
   1.483842,
   1.409255,
   1.476938,
   1.481399,
   1.482469,
   1.472752,
   1.485507,
   1.409255,
   1.409255,
   1.109894,
   1.481118,
   1.47787,
   1.488119,
   1.409255,
   1.480373,
   1.490046,
   1.481002,
   1.409255,
   1.505674,
   1.483522,
   1.48737,
   1.474276,
   1.483492,
   1.488044,
   1.497993,
   1.480387,
   1.485052,
   1.482208,
   1.487974,
   1.481675,
   1.492405,
   1.480498,
   1.48207,
   1.478592,
   1.479282,
   1.476754,
   1.490839,
   1.476254,
   1.472551,
   1.483077,
   1.477268,
   1.477091,
   1.486493,
   1.482491,
   1.479636,
   1.484816,
   1.483492,
   1.483492,
   1.508168,
   1.477445,
   1.437722,
   1.480222,
   1.489186,
   1.477891,
   1.491955,
   1.478859,
   1.484869,
   1.475427,
   1.48577,
   1.473861,
   1.489474,
   1.474881,
   1.480542,
   1.469344,
   1.475476,
   1.475314,
   1.485266,
   1.482516,
   1.480943,
   1.48017,
   1.473918,
   0.3025882,
   1.485093,
   1.48312,
   1.437722,
   1.475901,
   1.480293,
   1.437722,
   1.540275,
   1.514632,
   1.50951,
   1.503784,
   1.523824,
   1.514118,
   1.528196,
   1.503136,
   1.511589,
   1.503156,
   1.505787,
   1.500223,
   1.506599,
   1.503389,
   1.50794,
   1.498961,
   1.499368,
   1.497711,
   1.504144,
   1.494437,
   1.504345,
   1.500141,
   1.482031,
   1.503755,
   1.511857,
   1.518661,
   1.500911,
   1.504682,
   1.495952,
   1.504858,
   1.523548,
   1.497932,
   1.501839,
   1.512756,
   1.514672,
   1.512316,
   1.527073,
   1.5013,
   1.503089,
   1.50156,
   1.507681,
   1.498922,
   1.494209,
   1.500908,
   1.498477,
   1.496624,
   1.491059,
   1.503361,
   1.513841,
   1.503364,
   1.50322,
   1.504436,
   1.49998,
   1.502592,
   1.500888,
   1.514035,
   1.499954,
   1.505059,
   1.509297,
   1.507589,
   1.418912,
   1.418912,
   1.384417,
   1.481761,
   0.3755835,
   1.482817,
   1.498476,
   1.476687,
   1.472945,
   1.481101,
   1.418912,
   1.418912,
   1.488743,
   1.476877,
   1.418912,
   1.469282,
   1.475908,
   1.204811,
   1.486221,
   1.474391,
   1.480804,
   1.477968,
   1.480293,
   1.485306,
   1.485519,
   1.488016,
   1.477205,
   1.418912,
   1.471468,
   1.477282,
   1.495196,
   1.477589,
   1.492061,
   1.079365,
   1.135129,
   0.2098527,
   1.49581,
   1.474623,
   1.482216,
   1.192776,
   1.255245,
   1.192776,
   1.487961,
   1.478468,
   1.480454,
   1.47494,
   1.192776,
   0.2795169,
   0.1688917,
   1.487872,
   1.192776,
   1.245006,
   0.7230503,
   0.2581363,
   1.486025,
   1.490228,
   1.192776,
   1.480762,
   1.192776,
   1.488215,
   1.505206,
   1.447853,
   1.491149,
   1.470712,
   1.491457,
   1.482312,
   1.496168,
   1.479549,
   1.479473,
   1.447853,
   1.485806,
   1.479921,
   1.486614,
   1.485218,
   1.479943,
   1.447853,
   1.467752,
   1.475932,
   1.479224,
   1.48954,
   1.488547,
   1.490378,
   1.484943,
   1.482946,
   1.478363,
   0.4378605,
   1.499157,
   1.504858,
   1.497787,
   1.501208,
   1.489537,
   1.481116,
   1.480201,
   1.491495,
   1.13435,
   1.488985,
   1.487379,
   1.426824,
   1.467616,
   1.477729,
   1.480122,
   1.471462,
   1.482218,
   1.481236,
   1.481686,
   1.48209,
   1.4724,
   0.2739095,
   1.4871,
   1.47444,
   1.471898,
   1.477538,
   1.47552,
   1.47671,
   1.486533,
   1.479103,
   1.477862,
   1.479981,
   1.481722,
   1.485948,
   1.499612,
   1.477638,
   1.454786,
   1.487714,
   1.500904,
   1.492841,
   1.49564,
   1.469445,
   1.454786,
   1.129661,
   1.454786,
   1.454786,
   1.492352,
   1.479179,
   1.479057,
   1.479576,
   1.227704,
   1.454786,
   1.454786,
   1.454786,
   1.454786,
   1.485496,
   1.480271,
   1.481585,
   1.48975,
   1.484096,
   1.454786,
   1.478636,
   1.454786,
   1.484562,
   1.499193,
   1.483902,
   1.487216,
   1.482852,
   1.495722,
   1.483036,
   1.491829,
   1.476175,
   1.484487,
   1.476506,
   1.487075,
   1.486117,
   1.485665,
   1.477842,
   1.479101,
   1.479281,
   1.478675,
   1.47567,
   1.483036,
   1.478254,
   1.475023,
   1.479244,
   1.476234,
   1.485276,
   1.493422,
   1.482082,
   1.477317,
   1.482998,
   1.480229,
   1.487618,
   1.482012,
   1.488403,
   1.471013,
   1.489361,
   1.488261,
   1.476307,
   1.479803,
   1.487318,
   1.477356,
   1.482751,
   1.473913,
   1.468975,
   1.469738,
   1.475796,
   1.4675,
   1.469504,
   1.467113,
   1.477356,
   1.480357,
   1.482304,
   1.468366,
   1.46918,
   1.467998,
   1.477072,
   1.485779,
   1.484098,
   1.477356,
   1.484107,
   1.474242,
   1.477356,
   1.516061,
   1.506432,
   1.516048,
   1.507451,
   1.519037,
   1.511,
   1.514942,
   1.510572,
   1.498381,
   1.511181,
   1.512996,
   1.517024,
   0,
   0,
   0,
   0,
   0,
   0,
   1.50945,
   1.511307,
   1.496836,
   1.505455,
   1.509386,
   1.502104,
   1.498101,
   1.524039,
   1.502396,
   1.491246,
   1.503974,
   1.519219,
   1.524971,
   1.520291,
   1.51138,
   1.510075,
   1.510456,
   1.51806,
   1.516968,
   1.50062,
   1.506234,
   1.510123,
   1.508857,
   1.504062,
   0,
   0,
   0,
   0,
   0,
   0,
   1.509286,
   1.190005,
   1.505176,
   1.507596,
   1.507162,
   1.508796,
   1.52159,
   1.519148,
   1.514015,
   1.512564,
   1.504461,
   1.515086,
   1.495653,
   1.485791,
   1.479435,
   1.487465,
   1.485381,
   1.479435,
   1.479435,
   1.485801,
   1.484196,
   1.489824,
   1.484938,
   1.356582,
   0,
   0,
   0,
   0,
   0,
   0,
   1.486636,
   1.485906,
   1.479435,
   1.482628,
   1.479026,
   1.477613,
   1.488679,
   1.483521,
   1.483012,
   1.487536,
   1.48337,
   1.495139,
   1.496931,
   1.481145,
   1.396902,
   1.396902,
   1.396902,
   1.396902,
   1.495843,
   1.476801,
   1.484235,
   1.48219,
   1.396902,
   1.396902,
   1.487488,
   1.396902,
   0.3701772,
   1.470547,
   1.479867,
   1.478526,
   0.8970645,
   1.473378,
   1.396902,
   1.467441,
   1.396902,
   1.478261,
   1.396902,
   1.48593,
   1.477175,
   1.48682,
   1.480915,
   1.487312,
   1.520619,
   1.532687,
   1.548289,
   1.532726,
   1.544661,
   1.513736,
   1.497495,
   1.519077,
   1.513722,
   1.507187,
   1.516436,
   1.499597,
   1.499759,
   1.509567,
   1.496564,
   1.510174,
   1.48223,
   1.496908,
   1.51118,
   1.465341,
   1.491056,
   1.502914,
   1.497681,
   1.475426,
   1.517932,
   1.109235,
   1.482089,
   1.458423,
   1.491686,
   1.506182};
   TGraph *graph = new TGraph(540,tg_vdrift_vs_det_fx1,tg_vdrift_vs_det_fy1);
   graph->SetName("tg_vdrift_vs_det");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(0.6);
   
   TH1F *Graph_tg_vdrift_vs_det1 = new TH1F("Graph_tg_vdrift_vs_det1","",540,0,592.9);
   Graph_tg_vdrift_vs_det1->SetMinimum(0);
   Graph_tg_vdrift_vs_det1->SetMaximum(1.703118);
   Graph_tg_vdrift_vs_det1->SetDirectory(0);
   Graph_tg_vdrift_vs_det1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_tg_vdrift_vs_det1->SetLineColor(ci);
   Graph_tg_vdrift_vs_det1->GetXaxis()->SetLabelFont(42);
   Graph_tg_vdrift_vs_det1->GetXaxis()->SetLabelSize(0.035);
   Graph_tg_vdrift_vs_det1->GetXaxis()->SetTitleSize(0.035);
   Graph_tg_vdrift_vs_det1->GetXaxis()->SetTitleOffset(1);
   Graph_tg_vdrift_vs_det1->GetXaxis()->SetTitleFont(42);
   Graph_tg_vdrift_vs_det1->GetYaxis()->SetLabelFont(42);
   Graph_tg_vdrift_vs_det1->GetYaxis()->SetLabelSize(0.035);
   Graph_tg_vdrift_vs_det1->GetYaxis()->SetTitleSize(0.035);
   Graph_tg_vdrift_vs_det1->GetYaxis()->SetTitleFont(42);
   Graph_tg_vdrift_vs_det1->GetZaxis()->SetLabelFont(42);
   Graph_tg_vdrift_vs_det1->GetZaxis()->SetLabelSize(0.035);
   Graph_tg_vdrift_vs_det1->GetZaxis()->SetTitleSize(0.035);
   Graph_tg_vdrift_vs_det1->GetZaxis()->SetTitleOffset(1);
   Graph_tg_vdrift_vs_det1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_tg_vdrift_vs_det1);
   
   graph->Draw("p");
   
   Double_t tg_vdrift_vs_det_A_fx2[540] = {
   0,
   1,
   2,
   3,
   4,
   5,
   6,
   7,
   8,
   9,
   10,
   11,
   12,
   13,
   14,
   15,
   16,
   17,
   18,
   19,
   20,
   21,
   22,
   23,
   24,
   25,
   26,
   27,
   28,
   29,
   30,
   31,
   32,
   33,
   34,
   35,
   36,
   37,
   38,
   39,
   40,
   41,
   42,
   43,
   44,
   45,
   46,
   47,
   48,
   49,
   50,
   51,
   52,
   53,
   54,
   55,
   56,
   57,
   58,
   59,
   60,
   61,
   62,
   63,
   64,
   65,
   66,
   67,
   68,
   69,
   70,
   71,
   72,
   73,
   74,
   75,
   76,
   77,
   78,
   79,
   80,
   81,
   82,
   83,
   84,
   85,
   86,
   87,
   88,
   89,
   90,
   91,
   92,
   93,
   94,
   95,
   96,
   97,
   98,
   99,
   100,
   101,
   102,
   103,
   104,
   105,
   106,
   107,
   108,
   109,
   110,
   111,
   112,
   113,
   114,
   115,
   116,
   117,
   118,
   119,
   120,
   121,
   122,
   123,
   124,
   125,
   126,
   127,
   128,
   129,
   130,
   131,
   132,
   133,
   134,
   135,
   136,
   137,
   138,
   139,
   140,
   141,
   142,
   143,
   144,
   145,
   146,
   147,
   148,
   149,
   150,
   151,
   152,
   153,
   154,
   155,
   156,
   157,
   158,
   159,
   160,
   161,
   162,
   163,
   164,
   165,
   166,
   167,
   168,
   169,
   170,
   171,
   172,
   173,
   174,
   175,
   176,
   177,
   178,
   179,
   180,
   181,
   182,
   183,
   184,
   185,
   186,
   187,
   188,
   189,
   190,
   191,
   192,
   193,
   194,
   195,
   196,
   197,
   198,
   199,
   200,
   201,
   202,
   203,
   204,
   205,
   206,
   207,
   208,
   209,
   210,
   211,
   212,
   213,
   214,
   215,
   216,
   217,
   218,
   219,
   220,
   221,
   222,
   223,
   224,
   225,
   226,
   227,
   228,
   229,
   230,
   231,
   232,
   233,
   234,
   235,
   236,
   237,
   238,
   239,
   240,
   241,
   242,
   243,
   244,
   245,
   246,
   247,
   248,
   249,
   250,
   251,
   252,
   253,
   254,
   255,
   256,
   257,
   258,
   259,
   260,
   261,
   262,
   263,
   264,
   265,
   266,
   267,
   268,
   269,
   270,
   271,
   272,
   273,
   274,
   275,
   276,
   277,
   278,
   279,
   280,
   281,
   282,
   283,
   284,
   285,
   286,
   287,
   288,
   289,
   290,
   291,
   292,
   293,
   294,
   295,
   296,
   297,
   298,
   299,
   300,
   301,
   302,
   303,
   304,
   305,
   306,
   307,
   308,
   309,
   310,
   311,
   312,
   313,
   314,
   315,
   316,
   317,
   318,
   319,
   320,
   321,
   322,
   323,
   324,
   325,
   326,
   327,
   328,
   329,
   330,
   331,
   332,
   333,
   334,
   335,
   336,
   337,
   338,
   339,
   340,
   341,
   342,
   343,
   344,
   345,
   346,
   347,
   348,
   349,
   350,
   351,
   352,
   353,
   354,
   355,
   356,
   357,
   358,
   359,
   360,
   361,
   362,
   363,
   364,
   365,
   366,
   367,
   368,
   369,
   370,
   371,
   372,
   373,
   374,
   375,
   376,
   377,
   378,
   379,
   380,
   381,
   382,
   383,
   384,
   385,
   386,
   387,
   388,
   389,
   390,
   391,
   392,
   393,
   394,
   395,
   396,
   397,
   398,
   399,
   400,
   401,
   0,
   0,
   0,
   0,
   0,
   0,
   408,
   409,
   410,
   411,
   412,
   413,
   414,
   415,
   416,
   417,
   418,
   419,
   420,
   421,
   422,
   423,
   424,
   425,
   426,
   427,
   428,
   429,
   430,
   431,
   0,
   0,
   0,
   0,
   0,
   0,
   438,
   439,
   440,
   441,
   442,
   443,
   444,
   445,
   446,
   447,
   448,
   449,
   450,
   451,
   452,
   453,
   454,
   455,
   456,
   457,
   458,
   459,
   460,
   461,
   0,
   0,
   0,
   0,
   0,
   0,
   468,
   469,
   470,
   471,
   472,
   473,
   474,
   475,
   476,
   477,
   478,
   479,
   480,
   481,
   482,
   483,
   484,
   485,
   486,
   487,
   488,
   489,
   490,
   491,
   492,
   493,
   494,
   495,
   496,
   497,
   498,
   499,
   500,
   501,
   502,
   503,
   504,
   505,
   506,
   507,
   508,
   509,
   510,
   511,
   512,
   513,
   514,
   515,
   516,
   517,
   518,
   519,
   520,
   521,
   522,
   523,
   524,
   525,
   526,
   527,
   528,
   529,
   530,
   531,
   532,
   533,
   534,
   535,
   536,
   537,
   538,
   539};
   Double_t tg_vdrift_vs_det_A_fy2[540] = {
   1.472153,
   1.476872,
   1.329263,
   1.477937,
   1.46696,
   0.4685939,
   1.468299,
   1.471958,
   1.130933,
   1.466784,
   1.468952,
   1.473339,
   1.122207,
   1.469576,
   1.442505,
   1.189769,
   1.467284,
   1.329263,
   1.471428,
   1.466799,
   1.459458,
   1.473417,
   1.37242,
   1.483655,
   1.475213,
   1.337105,
   0.4844079,
   1.329263,
   1.3619,
   0.9701737,
   0.4206076,
   1.407309,
   1.407309,
   1.489943,
   1.477981,
   1.483859,
   1.407309,
   1.493605,
   1.462352,
   1.366873,
   1.407309,
   1.482901,
   1.483494,
   1.407309,
   1.467652,
   1.485096,
   1.477459,
   1.476519,
   1.485825,
   1.407309,
   1.407309,
   1.106048,
   1.474847,
   1.479352,
   1.490365,
   1.407309,
   1.4776,
   1.491859,
   1.479244,
   1.407309,
   1.477476,
   1.483997,
   1.462314,
   1.483099,
   1.4779,
   1.491864,
   1.476208,
   1.481611,
   1.463998,
   1.483588,
   1.467969,
   1.48158,
   1.478403,
   1.483353,
   1.466075,
   1.481169,
   1.472935,
   1.479103,
   1.479354,
   1.476383,
   1.465932,
   1.48444,
   1.474128,
   1.478725,
   1.480566,
   1.480488,
   1.473606,
   1.483646,
   1.4779,
   1.48919,
   1.475963,
   1.479197,
   1.432784,
   1.482893,
   1.463038,
   1.479113,
   1.470755,
   1.478416,
   1.463109,
   1.478341,
   1.466681,
   1.477579,
   1.475672,
   1.473767,
   1.465681,
   1.472245,
   1.466085,
   1.479769,
   1.474948,
   1.484352,
   1.471892,
   1.48154,
   1.471033,
   0.3312742,
   1.478621,
   1.483858,
   1.432784,
   1.477571,
   1.481787,
   1.432784,
   1.510474,
   1.516823,
   1.482285,
   1.506159,
   1.496336,
   1.518334,
   1.505615,
   1.503891,
   1.491736,
   1.50221,
   1.488362,
   1.502142,
   1.500766,
   1.503,
   1.499273,
   1.50109,
   1.490862,
   1.499771,
   1.494544,
   1.497869,
   1.498503,
   1.501223,
   1.477119,
   1.508269,
   1.505707,
   1.515688,
   1.498,
   1.503862,
   1.494602,
   1.508453,
   1.49658,
   1.49901,
   1.475458,
   1.515179,
   1.488833,
   1.514933,
   1.506146,
   1.502148,
   1.484094,
   1.502363,
   1.489225,
   1.500542,
   1.484668,
   1.500105,
   1.485314,
   1.496696,
   1.481823,
   1.501342,
   1.505563,
   1.502244,
   1.497234,
   1.504337,
   1.494615,
   1.503352,
   1.494011,
   1.514184,
   1.497344,
   1.502996,
   1.509283,
   1.50925,
   1.416363,
   1.416363,
   1.361359,
   1.487006,
   0.3746648,
   1.490041,
   1.477185,
   1.478556,
   1.454634,
   1.482583,
   1.416363,
   1.416363,
   1.47916,
   1.479534,
   1.416363,
   1.47218,
   1.473542,
   1.203445,
   1.479776,
   1.475013,
   1.480654,
   1.479654,
   1.481334,
   1.488219,
   1.478152,
   1.48735,
   1.473075,
   1.416363,
   1.474491,
   1.481094,
   1.46599,
   1.479955,
   1.464555,
   1.068944,
   1.113679,
   0.2173742,
   1.474293,
   1.476618,
   1.463003,
   1.185842,
   1.246502,
   1.185842,
   1.476746,
   1.479741,
   1.466516,
   1.478364,
   1.185842,
   0.280626,
   0.1737562,
   1.489136,
   1.185842,
   1.186455,
   0.7387578,
   0.272743,
   1.478466,
   1.490105,
   1.185842,
   1.484501,
   1.185842,
   1.493378,
   1.47902,
   1.442717,
   1.467974,
   1.473232,
   1.464248,
   1.48391,
   1.47767,
   1.481922,
   1.468626,
   1.442717,
   1.465608,
   1.481686,
   1.473135,
   1.487078,
   1.476695,
   1.442717,
   1.456738,
   1.481812,
   1.472798,
   1.488341,
   1.465557,
   1.494753,
   1.481598,
   1.482683,
   1.474186,
   0.4633164,
   1.500588,
   1.509489,
   1.494509,
   1.506198,
   1.492063,
   1.454165,
   1.484219,
   1.470216,
   1.141776,
   1.467824,
   1.487936,
   1.425671,
   1.47162,
   1.465203,
   1.483659,
   1.462456,
   1.486354,
   1.471723,
   1.481791,
   1.477775,
   1.477564,
   0.2959391,
   1.489432,
   1.46999,
   1.473235,
   1.480774,
   1.4801,
   1.478168,
   1.487337,
   1.478766,
   1.480122,
   1.483084,
   1.483153,
   1.488022,
   1.469516,
   1.481154,
   1.450909,
   1.491865,
   1.475945,
   1.495239,
   1.475771,
   1.472017,
   1.450909,
   1.13803,
   1.450909,
   1.450909,
   1.480742,
   1.480384,
   1.464428,
   1.483589,
   1.228293,
   1.450909,
   1.450909,
   1.450909,
   1.450909,
   1.485767,
   1.479243,
   1.48375,
   1.481971,
   1.482714,
   1.450909,
   1.480263,
   1.450909,
   1.487491,
   1.471879,
   1.487267,
   1.461367,
   1.489948,
   1.469082,
   1.477689,
   1.472982,
   1.476492,
   1.467654,
   1.478648,
   1.46876,
   1.48808,
   1.476772,
   1.479339,
   1.467505,
   1.479633,
   1.468677,
   1.479038,
   1.477689,
   1.482027,
   1.468921,
   1.481651,
   1.474603,
   1.48679,
   1.489309,
   1.481439,
   1.472858,
   1.4844,
   1.478721,
   1.491444,
   1.48415,
   1.463031,
   1.476665,
   1.464965,
   1.492127,
   1.457901,
   1.481681,
   1.465013,
   1.47115,
   1.464271,
   1.473695,
   1.457208,
   1.471478,
   1.460949,
   1.468887,
   1.457849,
   1.469458,
   1.47115,
   1.47852,
   1.471997,
   1.469539,
   1.460831,
   1.468319,
   1.476396,
   1.481448,
   1.47931,
   1.47115,
   1.479177,
   1.475026,
   1.47115,
   1.51628,
   1.482673,
   1.518476,
   1.479156,
   1.52177,
   1.48393,
   1.510448,
   1.489348,
   1.497105,
   1.490006,
   1.513647,
   1.500544,
   0,
   0,
   0,
   0,
   0,
   0,
   1.505151,
   1.495592,
   1.494682,
   1.491533,
   1.507341,
   1.494159,
   1.490107,
   1.512274,
   1.493162,
   1.484205,
   1.500885,
   1.515647,
   1.524739,
   1.489942,
   1.514138,
   1.483837,
   1.513936,
   1.490455,
   1.514395,
   1.479606,
   1.50628,
   1.49128,
   1.510445,
   1.485519,
   0,
   0,
   0,
   0,
   0,
   0,
   1.506027,
   1.190249,
   1.50492,
   1.49183,
   1.506703,
   1.497925,
   1.513369,
   1.50736,
   1.509201,
   1.501933,
   1.500136,
   1.509145,
   1.499483,
   1.459233,
   1.471658,
   1.463726,
   1.489901,
   1.471658,
   1.471658,
   1.467655,
   1.482001,
   1.470511,
   1.48762,
   1.337887,
   0,
   0,
   0,
   0,
   0,
   0,
   1.484965,
   1.474196,
   1.471658,
   1.47076,
   1.480005,
   1.474049,
   1.48293,
   1.475635,
   1.478866,
   1.480409,
   1.482253,
   1.49107,
   1.469797,
   1.4827,
   1.393459,
   1.393459,
   1.393459,
   1.393459,
   1.474013,
   1.478355,
   1.462962,
   1.484057,
   1.393459,
   1.393459,
   1.470817,
   1.393459,
   0.3873046,
   1.471135,
   1.470198,
   1.482359,
   0.904236,
   1.47158,
   1.393459,
   1.467679,
   1.393459,
   1.481469,
   1.393459,
   1.484411,
   1.471873,
   1.484695,
   1.480173,
   1.489372,
   1.494737,
   1.533185,
   1.52281,
   1.535777,
   1.518507,
   1.513778,
   1.478914,
   1.522619,
   1.499008,
   1.50704,
   1.501445,
   1.497401,
   1.493253,
   1.51203,
   1.4942,
   1.511169,
   1.476441,
   1.499449,
   1.506768,
   1.468106,
   1.488494,
   1.503609,
   1.492525,
   1.472955,
   1.513864,
   1.117916,
   1.479576,
   1.460863,
   1.487081,
   1.50892};
   graph = new TGraph(540,tg_vdrift_vs_det_A_fx2,tg_vdrift_vs_det_A_fy2);
   graph->SetName("tg_vdrift_vs_det_A");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(0.6);
   
   TH1F *Graph_tg_vdrift_vs_det_A2 = new TH1F("Graph_tg_vdrift_vs_det_A2","",540,0,592.9);
   Graph_tg_vdrift_vs_det_A2->SetMinimum(0);
   Graph_tg_vdrift_vs_det_A2->SetMaximum(1.689355);
   Graph_tg_vdrift_vs_det_A2->SetDirectory(0);
   Graph_tg_vdrift_vs_det_A2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_tg_vdrift_vs_det_A2->SetLineColor(ci);
   Graph_tg_vdrift_vs_det_A2->GetXaxis()->SetLabelFont(42);
   Graph_tg_vdrift_vs_det_A2->GetXaxis()->SetLabelSize(0.035);
   Graph_tg_vdrift_vs_det_A2->GetXaxis()->SetTitleSize(0.035);
   Graph_tg_vdrift_vs_det_A2->GetXaxis()->SetTitleOffset(1);
   Graph_tg_vdrift_vs_det_A2->GetXaxis()->SetTitleFont(42);
   Graph_tg_vdrift_vs_det_A2->GetYaxis()->SetLabelFont(42);
   Graph_tg_vdrift_vs_det_A2->GetYaxis()->SetLabelSize(0.035);
   Graph_tg_vdrift_vs_det_A2->GetYaxis()->SetTitleSize(0.035);
   Graph_tg_vdrift_vs_det_A2->GetYaxis()->SetTitleFont(42);
   Graph_tg_vdrift_vs_det_A2->GetZaxis()->SetLabelFont(42);
   Graph_tg_vdrift_vs_det_A2->GetZaxis()->SetLabelSize(0.035);
   Graph_tg_vdrift_vs_det_A2->GetZaxis()->SetTitleSize(0.035);
   Graph_tg_vdrift_vs_det_A2->GetZaxis()->SetTitleOffset(1);
   Graph_tg_vdrift_vs_det_A2->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_tg_vdrift_vs_det_A2);
   
   graph->Draw(" p");
   
   Double_t _fx3[394] = {
   0,
   0,
   1,
   3,
   4,
   5,
   6,
   7,
   8,
   9,
   10,
   11,
   12,
   13,
   14,
   16,
   18,
   19,
   20,
   21,
   22,
   23,
   24,
   25,
   26,
   28,
   29,
   42,
   44,
   45,
   46,
   47,
   60,
   61,
   62,
   63,
   65,
   66,
   67,
   68,
   69,
   70,
   71,
   72,
   73,
   74,
   75,
   76,
   77,
   78,
   79,
   80,
   81,
   82,
   83,
   84,
   85,
   86,
   87,
   89,
   90,
   91,
   93,
   94,
   95,
   96,
   97,
   98,
   99,
   100,
   101,
   102,
   103,
   104,
   105,
   106,
   107,
   108,
   109,
   110,
   111,
   112,
   120,
   121,
   122,
   123,
   124,
   125,
   126,
   127,
   128,
   129,
   130,
   131,
   133,
   134,
   135,
   136,
   137,
   138,
   139,
   140,
   141,
   142,
   143,
   144,
   145,
   146,
   147,
   148,
   149,
   150,
   151,
   152,
   153,
   154,
   155,
   156,
   157,
   158,
   159,
   160,
   161,
   162,
   163,
   164,
   165,
   166,
   167,
   168,
   169,
   170,
   171,
   172,
   173,
   174,
   175,
   176,
   177,
   178,
   179,
   186,
   187,
   188,
   189,
   192,
   193,
   195,
   196,
   197,
   198,
   199,
   200,
   201,
   202,
   203,
   204,
   205,
   206,
   208,
   209,
   210,
   211,
   212,
   214,
   216,
   217,
   218,
   220,
   222,
   223,
   224,
   225,
   229,
   232,
   240,
   242,
   243,
   244,
   245,
   246,
   247,
   248,
   250,
   251,
   252,
   253,
   254,
   256,
   257,
   258,
   259,
   260,
   261,
   262,
   263,
   264,
   265,
   266,
   267,
   268,
   269,
   270,
   271,
   272,
   273,
   274,
   275,
   276,
   278,
   279,
   280,
   281,
   282,
   283,
   284,
   285,
   286,
   288,
   289,
   290,
   291,
   292,
   293,
   294,
   295,
   296,
   297,
   298,
   299,
   300,
   301,
   303,
   304,
   305,
   312,
   313,
   314,
   315,
   316,
   330,
   331,
   332,
   333,
   334,
   336,
   337,
   338,
   339,
   340,
   341,
   342,
   343,
   344,
   345,
   346,
   347,
   349,
   350,
   351,
   352,
   353,
   354,
   355,
   356,
   357,
   358,
   359,
   360,
   361,
   362,
   363,
   364,
   365,
   366,
   367,
   369,
   370,
   371,
   372,
   373,
   374,
   375,
   376,
   378,
   379,
   380,
   381,
   382,
   383,
   390,
   391,
   392,
   393,
   394,
   395,
   396,
   397,
   398,
   399,
   400,
   401,
   408,
   409,
   410,
   411,
   412,
   413,
   414,
   415,
   416,
   417,
   418,
   419,
   420,
   421,
   422,
   423,
   424,
   425,
   426,
   427,
   428,
   429,
   430,
   431,
   438,
   439,
   440,
   441,
   442,
   443,
   444,
   445,
   446,
   447,
   448,
   449,
   457,
   458,
   459,
   460,
   468,
   469,
   471,
   472,
   473,
   474,
   475,
   476,
   477,
   478,
   479,
   486,
   487,
   488,
   489,
   492,
   494,
   495,
   496,
   497,
   505,
   506,
   507,
   508,
   509,
   510,
   511,
   512,
   513,
   514,
   515,
   516,
   517,
   518,
   519,
   520,
   521,
   522,
   523,
   524,
   525,
   526,
   527,
   528,
   529,
   530,
   531,
   532,
   534,
   535,
   536,
   539};
   Double_t _fy3[394] = {
   0,
   1.125742,
   1.135387,
   1.1099,
   1.050238,
   0.2842601,
   1.127975,
   1.19238,
   0.8545086,
   1.103227,
   1.03843,
   1.079242,
   0.7637146,
   1.158289,
   0.9655912,
   0.9230722,
   1.093142,
   1.027777,
   1.091538,
   1.036163,
   1.012437,
   1.093465,
   1.172522,
   1.023829,
   0.313804,
   0.9787604,
   0.520699,
   1.445346,
   1.378677,
   1.115835,
   1.102987,
   0.9835641,
   1.097643,
   1.036638,
   1.080532,
   1.139351,
   1.022898,
   1.014388,
   1.106164,
   1.094745,
   1.039248,
   1.007898,
   1.116809,
   0.8953306,
   1.050923,
   1.088615,
   1.115942,
   1.05222,
   1.118315,
   0.9681162,
   1.071512,
   1.094623,
   1.090828,
   1.053967,
   1.109407,
   0.9879754,
   1.054208,
   1.1375,
   1.194952,
   1.031209,
   1.004557,
   1.138215,
   1.132469,
   1.11754,
   1.054875,
   1.044642,
   1.084978,
   1.095216,
   1.192468,
   1.014719,
   1.072388,
   1.030908,
   1.077445,
   1.114784,
   1.066343,
   1.132046,
   1.135358,
   1.080378,
   1.078042,
   1.157575,
   1.11084,
   1.053264,
   1.390581,
   1.594771,
   1.567911,
   1.855941,
   1.03891,
   1.354373,
   0.91467,
   1.687036,
   1.382212,
   1.014051,
   1.572261,
   1.250825,
   1.322144,
   1.676087,
   1.6656,
   1.108321,
   1.411999,
   1.096365,
   1.335079,
   1.67653,
   1.575458,
   1.147013,
   1.226033,
   0.8492881,
   1.104028,
   1.484252,
   1.567873,
   1.172004,
   1.385814,
   1.118406,
   1.147778,
   0.9760624,
   1.002691,
   1.138713,
   1.066547,
   0.9981219,
   1.083381,
   1.059689,
   1.107428,
   1.075075,
   1.123778,
   0.9320637,
   1.046123,
   1.103455,
   1.146651,
   1.041689,
   1.091458,
   1.021984,
   1.036178,
   1.132306,
   1.172366,
   0.9784042,
   1.102388,
   0.6407317,
   1.119279,
   1.060842,
   1.15159,
   1.151961,
   1.18516,
   1.476214,
   1.565223,
   1.207873,
   1.195201,
   1.079746,
   1.052572,
   1.067334,
   1.066965,
   0.8120417,
   0.990221,
   0.9910718,
   1.189052,
   1.134631,
   1.016873,
   1.054556,
   0.8656214,
   1.200673,
   1.112529,
   1.066133,
   1.02688,
   1.038537,
   1.149974,
   1.130438,
   0.705751,
   1.077261,
   1.077836,
   0.9045608,
   1.344849,
   1.020639,
   1.159058,
   1.122075,
   1.130861,
   1.672603,
   0.8908015,
   1.413353,
   1.258298,
   1.211511,
   0.9253042,
   0.7163798,
   1.043704,
   1.01007,
   1.113646,
   1.020463,
   0.9519258,
   1.042327,
   1.090231,
   1.095457,
   0.9985347,
   0.8850498,
   0.9928582,
   1.106775,
   1.09858,
   1.140685,
   1.024118,
   0.935515,
   1.035762,
   0.2514199,
   1.119788,
   1.242683,
   1.233056,
   1.275268,
   1.026286,
   1.070999,
   1.094619,
   1.1431,
   0.8824473,
   0.9125266,
   1.087386,
   1.102986,
   1.00638,
   1.064676,
   0.9754016,
   1.055521,
   0.9763686,
   1.144426,
   1.034472,
   1.09014,
   1.080426,
   1.089027,
   1.060887,
   1.179787,
   1.064002,
   0.9603385,
   1.005943,
   1.02197,
   1.133917,
   1.089782,
   1.043327,
   1.075458,
   0.9913075,
   1.240733,
   1.06147,
   1.111181,
   1.148223,
   0.989636,
   1.089047,
   1.081553,
   1.084407,
   0.7803857,
   1.236541,
   1.132934,
   1.048359,
   1.049825,
   1.066832,
   1.018728,
   1.089664,
   1.055299,
   1.134146,
   1.145978,
   1.020469,
   0.9448333,
   1.043911,
   1.115779,
   1.096637,
   1.056435,
   1.044439,
   1.031305,
   1.049663,
   1.125064,
   1.086318,
   1.159573,
   0.937084,
   1.049295,
   1.09655,
   1.035163,
   0.9267398,
   0.9725907,
   1.012994,
   1.109984,
   1.046831,
   1.127018,
   1.036603,
   1.037709,
   0.9892314,
   1.305779,
   1.098208,
   1.007641,
   1.156262,
   1.092015,
   1.065387,
   1.123459,
   1.028923,
   1.02107,
   1.028606,
   1.107241,
   1.116049,
   1.126228,
   1.026587,
   0.9934285,
   0.9452108,
   1.043355,
   1.150121,
   1.508493,
   1.108692,
   1.350377,
   1.065139,
   1.623777,
   1.570295,
   1.493221,
   1.286828,
   1.033038,
   1.063394,
   1.598298,
   1.687317,
   1.591302,
   1.321452,
   1.103931,
   1.54266,
   1.344029,
   1.659855,
   1.346337,
   1.090873,
   1.240043,
   1.074275,
   1.084609,
   1.091391,
   1.152344,
   1.186619,
   1.130106,
   1.029072,
   1.025593,
   1.149246,
   1.18499,
   1.09127,
   1.030938,
   0.9276174,
   0.5590453,
   1.167655,
   1.089663,
   1.122007,
   1.197604,
   1.07045,
   1.020421,
   1.102231,
   1.154037,
   1.134791,
   0.9333523,
   1.008025,
   1.063669,
   1.16347,
   1.090438,
   1.19911,
   1.063994,
   1.092654,
   1.133279,
   1.052051,
   1.080921,
   1.095436,
   1.108313,
   1.076907,
   1.27176,
   1.000488,
   1.535858,
   1.493504,
   1.503142,
   1.27395,
   0.9046506,
   0.253604,
   1.124682,
   1.185134,
   1.116974,
   0.917302,
   1.17432,
   1.077062,
   1.040574,
   1.007152,
   1.008424,
   1.129427,
   1.218886,
   1.192161,
   1.157251,
   1.168976,
   1.022589,
   1.125814,
   1.155468,
   1.144104,
   1.151043,
   1.190457,
   1.0517,
   1.124685,
   1.024371,
   1.087611,
   1.083651,
   1.080791,
   1.004727,
   0.923797,
   1.123834,
   1.184452,
   1.048491,
   0.9481552,
   0.7168409,
   1.107312,
   1.214214};
   graph = new TGraph(394,_fx3,_fy3);
   graph->SetName("");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","",394,0,592.9);
   Graph_Graph3->SetMinimum(0);
   Graph_Graph3->SetMaximum(2.041535);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3->SetLineColor(ci);
   Graph_Graph3->GetXaxis()->SetTitle("drift velocity from OCDB, cm/us");
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleOffset(1);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetTitle("Detector ID");
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleOffset(1);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph3);
   
   graph->Draw("");
   tg_vdrift_vs_det_can->Modified();
   tg_vdrift_vs_det_can->cd();
   tg_vdrift_vs_det_can->SetSelected(tg_vdrift_vs_det_can);
}
